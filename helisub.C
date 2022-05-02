#include "src/exp_model.h"
#include "src/matrix1d.h"
#include "src/matrix2d.h"
#include "src/filename.h"
#include "src/metadata_table.h"

#include <iostream>


int main(int argc, char ** argv)
{

	if (argc < 7) 
	{
		std::cout << 
			"Usage: " << argv[0] << " <in> <out> <apix> <hrise> <hturn> <vectorX> <vectorY> <vectorZ> <offset> <#take> <boxsize> [<csym>]" <<
			std::endl ;
		exit(255);
	}

	FileName fn_in(argv[1]);
	FileName fn_out(argv[2]);

	RFLOAT apix = atof(argv[3]);
	RFLOAT hrise = atof(argv[4])/apix;
	RFLOAT hturn = atof(argv[5]);
	RFLOAT vx = atof(argv[6])/apix;
	RFLOAT vy = atof(argv[7])/apix;
	RFLOAT vz = atof(argv[8])/apix;
	RFLOAT offset = atof(argv[9]); 

	
	int take = atoi(argv[10]);

	int extract_size = atoi(argv[11]);

	int csym = 1;

	if (argc == 13) csym = atoi(argv[12]);

	Experiment my_data;

	my_data.read(fn_in);

	MetaDataTable MDout;

	FileName fn_img;

	int count=0;

	init_progress_bar(my_data.numberOfOriginalParticles());
	
	for (long int ipart=0; ipart< my_data.numberOfOriginalParticles(); ipart++)
	{

		//std::cout << ipart << " : rot, tilt, psi, ox, oy" << std::endl;
		MetaDataContainer * d = my_data.MDimg.getObject(ipart);
		
		RFLOAT rot, tilt, psi;
		RFLOAT ox, oy;
		my_data.MDimg.getValue(EMDL_ORIENT_ROT, rot, ipart);
		my_data.MDimg.getValue(EMDL_ORIENT_TILT, tilt, ipart);
		my_data.MDimg.getValue(EMDL_ORIENT_PSI, psi, ipart);
		my_data.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, ox, ipart);
		my_data.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, oy, ipart);

		//std::cout << "OLD: " << rot << " , " << tilt << " , " << psi << " , " << ox << " , " << oy << std::endl;

		Image<RFLOAT> Ipart, Isub;

		FileName fn_in;
		my_data.MDimg.getValue(EMDL_IMAGE_NAME, fn_in, ipart);

		Ipart.read(fn_in);

		//std::cout << "Filename: " << fn_in << " size : " << Ipart.getSize() << std::endl;
		//Ipart.write("test.mrcs", -1, true, WRITE_APPEND);

		for (int itake=0; itake < take; itake++)
		{
			RFLOAT dz, drot;
			int take2 = itake - take/2;
			dz = hrise * (take2 + offset);
			drot = hturn * (take2 + offset);
			
			Matrix1D<RFLOAT> shift, rotation;
			shift.resize(3);
			shift(0)=shift(1)=0;
			shift(2)=dz;

			rotation.resize(3);
			rotation(0)=rotation(1)=0;
			rotation(2)=drot;
		
				
			RFLOAT newrot, newtilt, newpsi;
			//newrot=rot+drot;
			//newtilt=tilt;
			//newpsi=psi;
			//



			Matrix1D<RFLOAT> point(3), new_point(3);
			point(0)=vx;
			point(1)=vy;
			point(2)=vz;

			Matrix2D<RFLOAT> helixrot(3,3), orientrot(3,3), finalrot(3,3);
			Euler_angles2matrix(drot,0,0, helixrot);
			Euler_angles2matrix(rot,tilt,psi, orientrot);

			new_point = helixrot * point + shift;
			
			new_point = orientrot * new_point;

			RFLOAT new_ox, new_oy;

			new_ox = ox + new_point(0);
			new_oy = oy + new_point(1);

			new_ox = - new_ox;
			new_oy = - new_oy;

			//std::cout << (int) floor(new_ox+0.5)-40 +360<< "\t" << (int) floor(new_oy+0.5)-40+360 << "\t80\t80\t-3" << std::endl;


			finalrot = orientrot * helixrot;
			Euler_matrix2angles(finalrot, newrot, newtilt, newpsi);
			

			//std::cout << newrot << " , " << newtilt << " , " << newpsi << " , " << new_ox << " , " << new_oy << std::endl;

			//extract particles
			int inewox, inewoy;
			inewox = (int) floor(new_ox);
			inewoy = (int) floor(new_oy);

			long int x0, xF, y0, yF;
			x0 = inewox + XSIZE(Ipart())/2 + FIRST_XMIPP_INDEX(extract_size);
	                xF = inewox + XSIZE(Ipart())/2 + LAST_XMIPP_INDEX(extract_size);
        	        y0 = inewoy + YSIZE(Ipart())/2 + FIRST_XMIPP_INDEX(extract_size);
               		yF = inewoy + YSIZE(Ipart())/2 + LAST_XMIPP_INDEX(extract_size);

		//	std::cout << "x0,xF,y0,yF: " << x0 << " , " << xF << " , " << y0 << " , " << yF << std::endl;

			Ipart().window(Isub(), y0, x0, yF, xF);


			int bg_radius;

			bg_radius = (extract_size / 2 - 1) ;

			normalise(Isub, bg_radius, -1, -1, true);

		//	fn_img.compose(count+1, "subparticles.mrcs");
			fn_img.compose(count+1, fn_out.getBaseName()+".mrcs");

			//std::cout << "Filename: " << fn_img << " size : " << Isub.getSize() << std::endl;
			Isub.write(fn_img, count+1, true, WRITE_APPEND);

			

			new_ox -=inewox;
			new_oy -=inewoy;

			MetaDataContainer * d2 = new MetaDataContainer(d->table, d);

			MDout.addObject(d2);
			MDout.setValue(EMDL_ORIENT_ROT, newrot);
			MDout.setValue(EMDL_ORIENT_TILT, newtilt);
			MDout.setValue(EMDL_ORIENT_PSI, newpsi);
			MDout.setValue(EMDL_ORIENT_ORIGIN_X, -new_ox);
			MDout.setValue(EMDL_ORIENT_ORIGIN_Y, -new_oy);
			MDout.setValue(EMDL_IMAGE_NAME, fn_img);

			count++;


		}



		if (ipart % 30 == 0 ) progress_bar(ipart);




	}

	MDout.write(fn_out);

	progress_bar(my_data.numberOfOriginalParticles());



}
