//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LandscapeConnectivityDriver.cpp
//
// Extraction of metrics related to landscape connectivity...
//	Get Dinf dirs and areas
//	Get Dinf catchment
//
// work in progress
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Martin D. Hurst
// British Geological Survey
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include "./../LSDStatsTools.hpp"
#include "./../LSDRaster.hpp"
#include "./../LSDIndexRaster.hpp"
#include "./../LSDFlowInfo.hpp"
#include "./../LSDJunctionNetwork.hpp"
#include "./../LSDRasterSpectral.hpp"
#include "./../TNT/tnt.h"

int main (int nNumberofArgs,char *argv[])
{
	//Test for correct input arguments
	if (nNumberofArgs!=3)
	{
		cout << "FATAL ERROR: wrong number inputs. The program needs the folder and DTM name" << endl;
		exit(EXIT_SUCCESS);
	}

	char TempChar[32];
	
	string Folder = argv[2];
	string DEM_name = argv[1];
	string dem_ext = "_dtm";
	string Fill_Name = "_fill";
	string Slope_Name = "_slope";
  	string DEM_flt_extension = "flt";
  	string DEM_txt_extension = "asc";
	
	//slopes to theshold
	//from np.logspace(-2,0,
	float slopes[] = {0.001,  0.0014, 0.0021, 0.003, 0.0043, 0.0062,  0.009,  0.013,  0.018,  0.026, 0.038,  0.055,  0.078,  0.11, 0.16, 0.23, 0.34,  0.48,  0.7,  1.};
	///float slopes[] = {0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.8,1.0};
	///float slopes[] = {0.1,0.5};
	vector<float> Slope_Thresholds(slopes, slopes+sizeof(slopes) / sizeof(float) );
	
	//load rasters in to memory
	LSDRaster Zeta((Folder+DEM_name+dem_ext), DEM_flt_extension);
	LSDRaster ZetaFill((Folder+DEM_name+Fill_Name), DEM_flt_extension);
	LSDRaster Slope((Folder+DEM_name+Slope_Name),DEM_flt_extension);
	
	//create an empty LSDRaster object for memory freeing
	int nrows = 1;
	int ncols = 1;
	float xmin = 0;
	float ymin = 0;
	float cellsize = 1;
	float ndv = -9999;
	Array2D<float> emptydata(nrows,ncols,ndv);
	LSDRaster Empty(nrows, ncols, xmin, ymin, cellsize, ndv, emptydata);
	
	//fill dtm
	//float MinSlope = 0.0001;		//minimum slope for dem fill function
	//LSDRaster ZetaFill = Zeta.fill(MinSlope);
	//string Fill_Name = "_fill";
	//ZetaFill.write_raster((Folder+DEM_name+Fill_Name),DEM_flt_extension);

	//vector<string> BoundaryConditions(4);
	//BoundaryConditions[0] = "No"; BoundaryConditions[1] = "No"; BoundaryConditions[2] = "No"; BoundaryConditions[3] = "No";
	//LSDFlowInfo FlowInfo(BoundaryConditions,ZetaFill); 
	//LSDRaster FlowLength = FlowInfo.distance_from_outlet();

	//write watershed to file
	//string FlowLength_Name = "_flow_length";
		
	//string Slope_Name = "_slope";
	//LSDRaster Slope((Folder+DEM_name+Slope_Name),DEM_flt_extension);
	
	//get slope raster
	//output extensions
//	vector<string> OutputExtensions;
//	OutputExtensions.push_back("_elev");
//	OutputExtensions.push_back("_slope");
//	OutputExtensions.push_back("_aspect");
//	OutputExtensions.push_back("_curv");
//	OutputExtensions.push_back("_pl_curv");
//	OutputExtensions.push_back("_pr_curv");
//	OutputExtensions.push_back("_tan_curv");
//	OutputExtensions.push_back("_sp");
//	
//	// Get surface metrics and output
//	vector<LSDRaster> SurfaceRasters;
//	vector<int> SelectedRasters(8,0);
//	SelectedRasters[1] = 1;
//	float WindowRadiusSmall = 10.;
//	
//	SurfaceRasters = Zeta.calculate_polyfit_surface_metrics(WindowRadiusSmall, SelectedRasters);
//	LSDRaster RasterToWrite;
//	
//	//loop through and write to file
//	for(int i=0, max=SelectedRasters.size(); i<max; ++i) {
//		if (SelectedRasters[i]==1) {
//      	RasterToWrite = SurfaceRasters[i];
//      	string OutputFile = Folder+DEM_name+OutputExtensions[i];
//      	RasterToWrite.write_raster((OutputFile), DEM_flt_extension);
////	}  }
//		
//	//assign slope
	//LSDRaster Slope = SurfaceRasters[1];

//	//setup an empty vector of LSDRasters to clear
//	vector<LSDRaster>().swap(SurfaceRasters);
	
	//FIRST GET AREA STATS FOR WHOLE CATCHMENT
	//compute Dinf and extract watershed
	cout << "Running Dinf...";
	Array2D<float> D_inf_dir_array = ZetaFill.D_inf_FlowDir();
	cout << "got D_inf_dir_array...";
	LSDRaster Dir = ZetaFill.write_dinf_flowdir_to_LSDRaster(D_inf_dir_array);
	cout << "created Dir LSDRaster...";
	LSDRaster Area = ZetaFill.D_inf_units();
	cout << "created Area LSDRaster...";
	cout << " Done!" << endl;
	
	//write area to file
	//string Area_Name = "_D_inf_area";
	//Area.write_raster((Folder+DEM_name+Area_Name),DEM_flt_extension);
	
	//write dir to file
	//string Dir_Name = "_D_inf_dir";
	//Dir.write_raster((Folder+DEM_name+Dir_Name),DEM_flt_extension);
	
	cout << "Getting Pour point...";
	float Amax = Get_Maximum(Area.get_RasterData(), Area.get_NoDataValue()); 
	int Index = Get_Minimum_Index(ZetaFill.get_RasterData(), ZetaFill.get_NoDataValue());
	int PourRow = (int)(Index/ZetaFill.get_NCols());
	int PourCol = Index % ZetaFill.get_NCols();
	cout << "Done!" << endl;
	
	//LSDIndexRaster D_inf_watershed = ZetaFill.D_inf_watershed(Dir, PourRow, PourCol);
	//int Catchment_Count = Get_Value_Count(D_inf_watershed.get_RasterData(), D_inf_watershed.get_NoDataValue());
	//float Catchment_Size = Catchment_Count*D_inf_watershed.get_DataResolution()*D_inf_watershed.get_DataResolution();
	cout << "Amax = " << Amax << "\n";
	
	//create mask of stream network for A > Acrit
	float Ac = 100000.;
	float Sc;
	string Condition = ">";
	LSDIndexRaster AreaMask = Area.Create_Mask(Condition,Ac);
	//LSDIndexRaster AreaMask2 = AreaMask.PadIndex();
	
	//write mask to file
	//string Area_Mask_Name = "AreaMask";
	//AreaMask.write_raster((Folder+DEM_name+Area_Mask_Name),DEM_txt_extension);

	//open text file for writing effective area results
	ofstream WriteResults("Results.txt");
	WriteResults << "Sc Amax AmaxMask\n";
	WriteResults << "0 0 0\n";
	
	Array2D<int> Mask;
	Array2D<int> SlopeMask; 
	
	//loop through slope values, threshold the DEM, create mask and rerun analysis
	for (int a=0, N=Slope_Thresholds.size(); a<N; ++a)
	{
		//get new copy of mask
		Mask = AreaMask.get_RasterData().copy();
		
		//get slope threshold
		Sc = Slope_Thresholds[a];
		cout << "*** Slope Threshold is " << Sc << endl;
		SlopeMask = Slope.Create_Mask(Condition, Sc).get_RasterData().copy();
		
//		LSDIndexRaster RawSlopeMask(AreaMask.get_NRows(), AreaMask.get_NCols(), AreaMask.get_XMinimum(), AreaMask.get_YMinimum(),
//            	AreaMask.get_DataResolution(), AreaMask.get_NoDataValue(), SlopeMask);

		//write mask to file
		//sprintf(TempChar,"_rawslopemask_%.3f",Sc);
		//string RSMask_Name = TempChar;
		//RawSlopeMask.write_raster((Folder+DEM_name+RSMask_Name),DEM_txt_extension);

		//combine the masks
		for (int i=0, NR=AreaMask.get_NRows(); i<NR; ++i) {
			for (int j=0, NC=AreaMask.get_NCols(); j<NC; ++j) {
				if (SlopeMask[i][j] == 1) Mask[i][j] = 1;
		}	}
	
		LSDIndexRaster TheMask(AreaMask.get_NRows(), AreaMask.get_NCols(), AreaMask.get_XMinimum(), AreaMask.get_YMinimum(),
            	AreaMask.get_DataResolution(), AreaMask.get_NoDataValue(), Mask);
	
		//write mask to file
		//sprintf(TempChar,"_mask_%.3f",Sc);
		//string Mask_Name = TempChar;
		//TheMask.write_raster((Folder+DEM_name+Mask_Name),DEM_txt_extension);
		
		//extract DTM from mask
		//LSDRaster ZetaMasked = Zeta.ExtractByMask(TheMask);
		//sprintf(TempChar,"_Z_mask_%.3f",Sc);
		//string Zeta_Mask_Name = TempChar;
		//ZetaMasked.write_raster((Folder+DEM_name+Zeta_Mask_Name),DEM_txt_extension);
		
		//Fill
		LSDRaster ZetaMaskedFill = ZetaFill.ExtractByMask(TheMask);
		
		//Get Dinf Dir, Area and Watershed
		LSDRaster DirMasked = Dir.ExtractByMask(TheMask);
		//LSDRaster AreaMasked = Area.ExtractByMask(TheMask);
				
		//Array2D<float> D_inf_dir_array = ZetaMaskedFill.D_inf_FlowDir();
		//LSDRaster Dir = ZetaMaskedFill.write_dinf_flowdir_to_LSDRaster(D_inf_dir_array);
		LSDRaster FlowAcc = ZetaMaskedFill.D_inf_FlowArea(DirMasked.get_RasterData());
		LSDRaster AreaMasked = FlowAcc.D_inf_ConvertFlowToArea();

		//write area to file
		//sprintf(TempChar, "_Area_%.3f", Sc);
		//string Area_Name = TempChar;
		//AreaMasked.write_raster((Folder+DEM_name+Area_Name),DEM_flt_extension);
	
		//write dir to file
		//sprintf(TempChar, "_Dir_%.3f", Sc);
		//string Dir_Name = TempChar;
		//DirMasked.write_raster((Folder+DEM_name+Dir_Name),DEM_flt_extension);
	
		float Amax_Mask = Get_Maximum(AreaMasked.get_RasterData(), AreaMasked.get_NoDataValue()); 
		//int Index = Get_Minimum_Index(ZetaMaskedFill.get_RasterData(), ZetaMaskedFill.get_NoDataValue());
		//int PourRow = (int)(Index/ZetaFill.get_NCols());
		//int PourCol = Index % ZetaFill.get_NCols();
	
		LSDIndexRaster D_inf_watershed = ZetaMaskedFill.D_inf_watershed(DirMasked, PourRow, PourCol);
		
		sprintf(TempChar,"_Ae_%.3f",Sc);
		string Ae_Name = TempChar;
		D_inf_watershed.write_raster((Folder+DEM_name+Ae_Name),DEM_flt_extension);

		//int Catchment_Count = Get_Value_Count(D_inf_watershed.get_RasterData(), D_inf_watershed.get_NoDataValue());
		//float Catchment_Size_Mask = Catchment_Count*D_inf_watershed.get_DataResolution()*D_inf_watershed.get_DataResolution();
		
		WriteResults << Sc << " " << Amax << " " << Amax_Mask << "\n";
	}
	WriteResults.close();
}
