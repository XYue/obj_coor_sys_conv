#include <iostream>
#include <fstream>
#include <iomanip>

#include <Eigen/Dense>
#include <proj_api.h>
#include <tinyxml2.h>

using namespace tinyxml2;

int ReadFromXml(const std::string & param_file,
				Eigen::Vector3d & translation,
				Eigen::Matrix3d & rotation,
				double & scale,
				std::string & dproj_cmd)
{
	int rtn = -1;

	do 
	{
#undef  SAFE_TEXT
#define SAFE_TEXT(txt)  (txt ? txt : "")

		translation = Eigen::Vector3d::Zero();
		rotation = Eigen::Matrix3d::Zero();
		scale = 0.;
		dproj_cmd = "";

		Eigen::Vector3d t = Eigen::Vector3d::Zero();
		Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
		double s = 0.;
		std::string cmd = "";


		XMLDocument doc;
		if (doc.LoadFile(param_file.c_str())!= XML_NO_ERROR)
			goto error0;

		// search module list
		const XMLElement *module_list = doc.FirstChildElement("module_list");
		if(!module_list)
			goto error0;

		// search module item
		bool read_all = false;
		for(const XMLElement *module_item = module_list->FirstChildElement("module_item");
			module_item != NULL; module_item = module_item->NextSiblingElement("module_item"))
		{
			t = Eigen::Vector3d::Zero();
			R = Eigen::Matrix3d::Zero();
			s = 0.;
			cmd = "";

			// search module name & version
			const XMLElement *module_name = module_item->FirstChildElement("module_name");
			const XMLElement *module_version = module_item->FirstChildElement("module_version");
			
			if(std::string("mesh_painting").compare(SAFE_TEXT(module_name->GetText())))
			{
				continue;
			}
			if(std::string("1.0").compare(SAFE_TEXT(module_version->GetText())))
			{
				continue;
			}


			// 处理参数
			const XMLElement *mp_param = 
				module_item->FirstChildElement("model_projection_param");
			if (!mp_param) continue;
						
			const XMLElement * aom_param =
				mp_param->FirstChildElement("aom_parameters");
			if (aom_param)
			{
				std::stringstream sstr;

				// scale
				const XMLElement * scale = 
					aom_param->FirstChildElement("scale");
				if(scale)
				{
					double test = atof(SAFE_TEXT(scale->GetText()));
					if (test > 0.)
					{
						s = test;
					}
				} else continue;

				// translation
				const XMLElement *translation = 
					aom_param->FirstChildElement("translation");
				if(translation)
				{ 
					Eigen::Vector3d tmp_t;
					unsigned int c = 0;
					sstr.clear(); sstr.str("");
					sstr << SAFE_TEXT(translation->GetText());
					for(c = 0; c < 3 && sstr >> tmp_t(c); c++) {}
					if(c != 3)
					{
						continue;
					}
					t = tmp_t;
				}

				// rotation
				const XMLElement *rot_mat = 
					aom_param->FirstChildElement("rotation");
				if(rot_mat)
				{
					Eigen::Matrix3d tmp_R;
					unsigned int c = 0;
					unsigned int r = 0;
					for(const XMLElement *row = 
						rot_mat->FirstChildElement("row");
						r < 3 && row != NULL; 
					r++, row = row->NextSiblingElement("row"))
					{
						sstr.clear(); sstr.str("");
						sstr << SAFE_TEXT(row->GetText());
						for(c = 0; c < 3 && sstr >> tmp_R(r,c); c++) {}
						if(c < 3) break;								
					}
					if(r != 3 || c != 3)
					{
						continue;
					} 

					R = tmp_R;
				}

				const XMLElement *xml_dproj_cmd = 
					mp_param->FirstChildElement("dproj_cmd");
				if (xml_dproj_cmd)
				{
					std::string test = SAFE_TEXT(xml_dproj_cmd->GetText());
					if (!test.empty())
					{
						cmd = test;
					}
				}

				read_all = true;
			}
		}

		if (read_all)
		{
			translation = t;
			rotation = R;
			scale = s;
			dproj_cmd = cmd;
		} else goto error0;

		rtn = 0;
	} while (0);
error0:
	
	return rtn;
}

int utm_to_GK3(const std::string & proj_cmd,
			   const Eigen::Vector3d & utm,
			   Eigen::Vector3d & gk3,
			   bool prefix_zone_code = true)
{
	int rtn = -1;

	do 
	{
		gk3(1) = utm(1) / 0.9996;
		gk3(0) = (utm(0) - 500000.) / 0.9996 + 500000;
		gk3(2) = utm(2);

		if (prefix_zone_code && !proj_cmd.empty())
		{
			std::string dst_prj_cmd = "+proj=longlat +datum=WGS84 +ellps=WGS84";
			std::string src_prj_cmd = proj_cmd/*"+proj=utm +zone=49 +datum=WGS84 +ellps=WGS84"*/;
			projPJ sproj = pj_init_plus(src_prj_cmd.c_str());
			projPJ dproj = pj_init_plus(dst_prj_cmd.c_str());
			if(sproj && dproj) 
			{
				double longlat[3];
				longlat[0] = utm[0];
				longlat[1] = utm[1];
				longlat[2] = utm[2];

				if(pj_is_latlong(sproj))
				{
					goto error0;
				}

				if(!pj_transform(sproj, dproj, 1, 0, longlat, longlat+1, longlat+2))
				{
					if(pj_is_latlong(dproj))
					{
						longlat[0] *= RAD_TO_DEG;
						longlat[1] *= RAD_TO_DEG;
					} else goto error0;
				} else goto error0;

				int zone = floor(floor((longlat[0]-1.5)/3.)+1);
				int pa = 1;
				while (static_cast<int>(gk3[0] / pa) > 0)
				{
					pa *= 10;
				}
				gk3[0] += pa * zone;
				
			} else goto error0;
		}

		rtn = 0;
	} while (0);
error0:

	return rtn;
}

void main(int argc, char **argv)
{
	int rtn = -1;

	do
	{
		bool output_gk3 = false;

		if (argc != 4)
		{
			std::cout<<"Usage:"<<std::endl;
			std::cout<<"\tobj_coor_sys_conv.exe input_obj_filename param_input_xml_filename output_obj_filename"<<std::endl;
			std::cout<<"\te.g. obj_coor_sys_conv.exe mesh.obj input.xml out.obj"<<std::endl;
			break;
		}

		std::string input_file = argv[1];
		std::string param_file = argv[2];
		std::string output_file = argv[3];
		std::string line, word;
		std::stringstream sstr;

		Eigen::Matrix3d R;
		Eigen::Vector3d translation;
		double scale;
		std::string proj_cmd;
		
		if (ReadFromXml(param_file, translation, R, scale, proj_cmd))
		{
			std::cout<<"Loading parameters from "<<param_file<<" failed."<<std::endl;
			goto error0;
		}

		std::ifstream i_file(input_file);
		std::ofstream o_file(output_file);
		if (!i_file.good() || !o_file.good())
		{
			break;
		}

		o_file << std::fixed;
		o_file << std::setprecision(15);
		while (i_file.good() && o_file.good())
		{
			word.clear();word = "";
			line.clear();line = "";
			sstr.clear();sstr.str("");

			getline(i_file, line);
			sstr << line;
			sstr >> word;
			if (word.compare("v"))
			{
				//std::cout<<line<<std::endl;
				o_file << line << std::endl;
				continue;
			}

			Eigen::Vector3d ori_v;
			Eigen::Vector3d utm;
			Eigen::Vector3d dst_v;

			sstr >> ori_v(0) >> ori_v(1) >> ori_v (2);
			utm = scale * R * ori_v;

			if (output_gk3)
			{
				utm += translation;

				if (utm_to_GK3(proj_cmd, utm, dst_v,1))
				{
					std::cout<<"Converting from utm to GK3 failed."<<std::endl;
					goto error0;
				}
			} else {
				dst_v = utm;
			}

			o_file << "v " 
				<< dst_v(0) << " "
				<< dst_v(1) << " "
				<< dst_v(2) << std::endl;
		}

		i_file.close();
		o_file.close();


		// save translation point
		if (!output_gk3)
		{
			std::string longlat_proj_cmd = "+proj=longlat +datum=WGS84 +ellps=WGS84";
			projPJ sproj = pj_init_plus(proj_cmd.c_str());
			projPJ dproj = pj_init_plus(longlat_proj_cmd.c_str());
			bool valid = false;

			do 
			{
				if(!sproj || !dproj) break;

				double longlat[3];
				longlat[0] = translation[0];
				longlat[1] = translation[1];
				longlat[2] = translation[2];

				if(pj_is_latlong(sproj)) break;

				if(!pj_transform(sproj, dproj, 1, 0, longlat, longlat+1, longlat+2) && 
					pj_is_latlong(dproj))
				{
					longlat[0] *= RAD_TO_DEG;
					longlat[1] *= RAD_TO_DEG;
				} else break;

				std::ofstream longlat_file(output_file + ".longlat");
				if (longlat_file.good())
				{
					longlat_file << std::setprecision(15)
						<< longlat[0] << " " << longlat[1] << " " << longlat[2] << std::endl;
					longlat_file.close();
				}

				valid = true;
			} while (0);

			if(dproj)
			{
				pj_free(dproj);
				dproj = NULL;
			}
			if(sproj)
			{
				pj_free(sproj);
				sproj = NULL;
			}

			if (!valid) break;
		}
		


		rtn = 0;

	}while (0);
error0:

	if (rtn)
		std::cout<<"fail."<<std::endl;
	else 
		std::cout<<"done."<<std::endl;
}