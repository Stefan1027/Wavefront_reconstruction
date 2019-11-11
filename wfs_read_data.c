/*===============================================================================================================================

	Date:          Sep-20-18
	Software-Nr:   N/A
	Version:       1.3
	Author:        Stefan-Tudor Ilie (sti1g14@soton.ac.uk)

	Changelog:     Aug-23-18 -> 1.0 -> First version created
				   Sep-05-18 -> 1.1 -> Updated script to generalise for all WFS cameras and created default sensor settings
				   Sep-20-18 -> 1.2 -> Fixed delete row/column functions and added average set of measurements
				   Sep-20-18 -> 1.3 -> Corrected wavefront geometry issues, and outputed both geometries (Pixel and um dimensions)


===============================================================================================================================*/


/*===============================================================================================================================
  Include Files

  Note: You may need to set your compilers include search path to the VXIPNP include directory.
		  This is typically 'C:\Program Files (x86)\IVI Foundation\VISA\WinNT\WFS'.

===============================================================================================================================*/

#include "WFS.h" // Wavefront Sensor driver's header file
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <utility.h>
//#include <userint.h>



/*===============================================================================================================================
  Defines
===============================================================================================================================*/

#define  DEVICE_OFFSET_WFS10           (0x00100) // device IDs of WFS10 instruments start at 256 decimal
#define  DEVICE_OFFSET_WFS20           (0x00200) // device IDs of WFS20 instruments start at 512 decimal
#define  DEVICE_OFFSET_WFS30           (0x00400) // device IDs of WFS30 instruments start at 1024 decimal
#define  DEVICE_OFFSET_WFS40           (0x00800) // device IDs of WFS40 instruments start at 2048 decimal

// settings for this sample program, you may adapt settings to your preferences
#define  OPTION_OFF                    (0)
#define  OPTION_ON                     (1)

#define  SAMPLE_PIXEL_FORMAT           PIXEL_FORMAT_MONO8   // only 8 bit format is supported
#define  SAMPLE_REF_PLANE              WFS_REF_INTERNAL
#define  SAMPLE_IMAGE_READINGS         (10) // trials to read a exposed spotfield image

#define  SAMPLE_OPTION_DYN_NOISE_CUT   OPTION_ON   // use dynamic noise cut features  
#define  SAMPLE_OPTION_CALC_SPOT_DIAS  OPTION_OFF  //  calculate spot diameters
#define  SAMPLE_OPTION_CANCEL_TILT     OPTION_OFF   // do not cancel average wavefront tip and tilt
#define  SAMPLE_OPTION_LIMIT_TO_PUPIL  OPTION_OFF  // don't limit wavefront calculation to pupil interior

#define  SAMPLE_OPTION_HIGHSPEED       OPTION_ON   // use highspeed mode (only for WFS10 and WFS20 instruments)
#define  SAMPLE_OPTION_HS_ADAPT_CENTR  OPTION_ON   // adapt centroids in highspeed mode to previously measured centroids
#define  SAMPLE_HS_NOISE_LEVEL         (30)        // cut lower 30 digits in highspeed mode
#define  SAMPLE_HS_ALLOW_AUTOEXPOS     (1)         // allow autoexposure in highspeed mode (runs somewhat slower)

#define  SAMPLE_WAVEFRONT_TYPE         WAVEFRONT_MEAS // calculate measured wavefront

// Normal outputs, saved as reference
#define  OUTPUT_FILE_NAME        		 	"WFS_output_file.csv"
#define  CENTROID_DEVIATION_FILE_NAME_X  	"WFS_deviation_file_x.csv"
#define  CENTROID_DEVIATION_FILE_NAME_Y 	"WFS_deviation_file_y.csv"
#define  FIELD_INTENSITY_FILE_NAME   		"WFS_field_intensity.csv"
#define  REFERENCE_SPOT_POS_X   	     	"WFS_reference_spot_x.csv"
#define  REFERENCE_SPOT_POS_Y   	     	"WFS_reference_spot_y.csv"
#define  CENTROID_POS_X   	    		 	"WFS_centroid_pos_x.csv"
#define  CENTROID_POS_Y   	    		 	"WFS_centroid_pos_y.csv"
#define  CENTROID_DEVIATION_FILE_NAME_X_AVG "WFS_deviation_file_x_avg.csv"
#define  CENTROID_DEVIATION_FILE_NAME_Y_AVG "WFS_deviation_file_y_avg.csv"
#define  WAVEFRONT_ZERNIKE_FIT		     	"WFS_wavefront_zernike_fit.csv"
#define  FIELD_INTENSITY_FILE_AVG  			"WFS_field_intensity_avg.csv"
#define  MLA_FOCAL							"WFS_focal.csv"

// Processed outputs, based on selective row & column removal

#define  PROCESS_CENTROID_DEVIATION_FILE_NAME_X  		"WFS_deviation_file_x_processed.csv"
#define  PROCESS_CENTROID_DEVIATION_FILE_NAME_Y  		"WFS_deviation_file_y_processed.csv"
#define  PROCESS_FIELD_INTENSITY_FILE_NAME   	 		"WFS_field_intensity_processed.csv"
#define  PROCESS_REFERENCE_SPOT_POS_X   	     		"WFS_reference_spot_x_processed.csv"
#define  PROCESS_REFERENCE_SPOT_POS_Y   	     		"WFS_reference_spot_y_processed.csv"
#define  PROCESS_CENTROID_POS_X   	    		 		"WFS_centroid_pos_x_processed.csv"
#define  PROCESS_CENTROID_POS_Y   	    		 		"WFS_centroid_pos_y_processed.csv"
#define  PROCESS_CENTROID_DEVIATION_FILE_NAME_X_AVG   "WFS_deviation_file_x_avg_processed.csv"
#define  PROCESS_CENTROID_DEVIATION_FILE_NAME_Y_AVG   "WFS_deviation_file_y_avg_processed.csv"
#define  PROCESS_WAVEFRONT_ZERNIKE_FIT		     		"WFS_wavefront_zernike_fit_processed.csv"
#define  PROCESS_FIELD_INTENSITY_FILE_AVG  				"WFS_field_intensity_avg_processed.csv"



// max_spots = Resolution * pixel_size / lenslet pitch  
// lenslet pitch = 150 um by default for 5C and 7AR lenslets and 300um for 14AR versions
// 1080 * 5 /150 = 36 maximum spots at 1080x1080 resolution

#ifdef __cplusplus
    extern "C" {
#endif


/*===============================================================================================================================
  Data type definitions
===============================================================================================================================*/

typedef struct
{
	int               selected_id;
	int               handle;
	int               status;
	
	char              version_wfs_driver[WFS_BUFFER_SIZE];
	char              version_cam_driver[WFS_BUFFER_SIZE];
	char              manufacturer_name[WFS_BUFFER_SIZE];
	char              instrument_name[WFS_BUFFER_SIZE];
	char              serial_number_wfs[WFS_BUFFER_SIZE];
	char              serial_number_cam[WFS_BUFFER_SIZE];
	
	int               mla_cnt;
	int               selected_mla;
	int               selected_mla_idx;
	char              mla_name[WFS_BUFFER_SIZE];
	double            cam_pitch_um;
	double            lenslet_pitch_um;
	double            center_spot_offset_x;
	double            center_spot_offset_y;
	double            lenslet_f_um;
	double            grd_corr_0;
	double            grd_corr_45;
	int               spots_x;
	int               spots_y;
	

}  instr_t;


/*===============================================================================================================================
  Function Prototypes
===============================================================================================================================*/
void handle_errors (int);
int select_instrument (int *selection, ViChar resourceName[]);
int select_mla (int *selection);
void deleteRow(float array[], int row, int rows, int cols);
void deleteCol(float array[], int col, int rows, int cols);

/*===============================================================================================================================
  Global Variables
===============================================================================================================================*/

const int   zernike_modes[] = { 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66 }; // converts Zernike order to Zernike modes
instr_t     instr = { 0 };    // all instrument related data are stored in this structure

int         hs_win_count_x,hs_win_count_y,hs_win_size_x,hs_win_size_y; // highspeed windows data
int         hs_win_start_x[MAX_SPOTS_X],hs_win_start_y[MAX_SPOTS_Y];


/*===============================================================================================================================
  Code
===============================================================================================================================*/
void main ()
{
	int               err=0;
	int               i,j,cnt;
	int               rows, cols;   // image height and width, depending on camera resolution
	int               selection;
	int 			  iterations;
	unsigned char     *ImageBuffer; // pointer to the camera image buffer
	
	double            expos_act, master_gain_act;
	double            beam_centroid_x, beam_centroid_y;
	double            beam_diameter_x, beam_diameter_y;
	double		      pupil_centre_x,   pupil_centre_y;
	double			  pupil_diameter_x, pupil_diameter_y;
	float             centroid_x[MAX_SPOTS_Y][MAX_SPOTS_X];
	float             centroid_y[MAX_SPOTS_Y][MAX_SPOTS_X];
    float             intensity[MAX_SPOTS_X][MAX_SPOTS_Y];

	float             deviation_x[MAX_SPOTS_Y][MAX_SPOTS_X];
	float             deviation_y[MAX_SPOTS_Y][MAX_SPOTS_X];
	float 			  deviation_avg_x[MAX_SPOTS_Y][MAX_SPOTS_X];
	float 			  deviation_avg_y[MAX_SPOTS_Y][MAX_SPOTS_X];
	float             wavefront[MAX_SPOTS_Y][MAX_SPOTS_X];
	float			  intensity_avg[MAX_SPOTS_Y][MAX_SPOTS_X];
	
	float             zernike_um[MAX_ZERNIKE_MODES+1];             // index runs from 1 - MAX_ZERNIKE_MODES
	float             zernike_orders_rms_um[MAX_ZERNIKE_ORDERS+1]; // index runs from 1 - MAX_ZERNIKE_MODES
	float 	 		  spot_ref_x[MAX_SPOTS_Y][MAX_SPOTS_X];
	float 	 		  spot_ref_y[MAX_SPOTS_Y][MAX_SPOTS_X];
	double            roc_mm;
	double			  spot_min,spot_max,spot_mean;
	
	int               zernike_order;
	int 			  rows_columns_del;
	float			  SAMPLE_PUPIL_CENTROID_Y, SAMPLE_PUPIL_CENTROID_X, SAMPLE_PUPIL_DIAMETER_X, SAMPLE_PUPIL_DIAMETER_Y;
	int 			  SAMPLE_ZERNIKE_ORDERS, SAMPLE_PRINTOUT_SPOTS;
	int   			  device_resolution;
	int 			  ROWS,COLS;
	

	double            wavefront_min, wavefront_max, wavefront_diff, wavefront_mean, wavefront_rms, wavefront_weighted_rms;

	ViChar            resourceName[256];
	FILE              *fp;
    FILE              *fp1;
    FILE              *fp2;
    FILE              *fp3;
	FILE 		      *fp4;
	FILE		      *fp5;
	FILE			  *fp6;
	FILE			  *fp7;
	FILE			  *fp8;
	FILE              *fp9;
	FILE              *fp10;
	FILE 			  *fp11;
	FILE  			  *fp12;


	printf("This is a Thorlabs Wavefront Sensor application.\n\n");
	
	// Get the driver revision
	if(err == WFS_revision_query (NULL, instr.version_wfs_driver, instr.version_cam_driver)) // pass NULL because handle is not yet initialized
		handle_errors(err);
	
	//printf("Camera USB driver version     : %s\n", instr.version_cam_driver);
	printf("WFS instrument driver version : %s\n\n", instr.version_wfs_driver);
	
	
	// Show all and select one WFS instrument
	if(select_instrument(&instr.selected_id, resourceName) == 0)
	{
		printf("\nNo instrument selected. Press <ENTER> to exit.\n");
		fflush(stdin);
		getchar();
		return; // program ends here if no instrument selected
	}
	
	// Get the resource name for this instrument
	//if(err == WFS_GetInstrumentListInfo (VI_NULL, instr.selected_id, VI_NULL, VI_NULL, VI_NULL, VI_NULL, resourceName))
	// handle_errors(err);
	
	
	// print out the resource name
	printf("\nResource name of selected WFS: %s\n", resourceName);
	
	
	// Open the Wavefront Sensor instrument
	//if(err == WFS_init (instr.selected_id, &instr.handle))
	if(err == WFS_init (resourceName, VI_FALSE, VI_FALSE, &instr.handle)) 
		handle_errors(err);

	// Get instrument information
	if(err == WFS_GetInstrumentInfo (instr.handle, instr.manufacturer_name, instr.instrument_name, instr.serial_number_wfs, instr.serial_number_cam))
		handle_errors(err);
	
	printf("\n");
	printf("Opened Instrument:\n");
	printf("Manufacturer           : %s\n", instr.manufacturer_name);
	printf("Instrument Name        : %s\n", instr.instrument_name);
	printf("Serial Number WFS      : %s\n", instr.serial_number_wfs);
	
	
	// Select a microlens array (MLA)
	if(select_mla(&instr.selected_mla) < 0)
	{
		printf("\nNo MLA selected. Press <ENTER> to exit.\n");
		fflush(stdin);
		getchar();
		return;
	}
	
	// Activate desired MLA
	if(err == WFS_SelectMla (instr.handle, instr.selected_mla))
		handle_errors(err);



	printf("\nInsert pupil centroid x coordinate\n");
	fflush(stdin);
	scanf("%f", &SAMPLE_PUPIL_CENTROID_X);

	printf("\nInsert pupil centroid y coordinate\n");
	fflush(stdin);
	scanf("%f", &SAMPLE_PUPIL_CENTROID_Y);	

	printf("\nInsert pupil diameter x coordinate\n");
	fflush(stdin);
	scanf("%f", &SAMPLE_PUPIL_DIAMETER_X);

	printf("\nInsert pupil diameter y coordinate\n");
	fflush(stdin);
	scanf("%f", &SAMPLE_PUPIL_DIAMETER_Y);

	printf("\nInsert number of zernike orders diameter to be fitted to the wavefront (max 11)\n");
	fflush(stdin);
	scanf("%d", &SAMPLE_ZERNIKE_ORDERS);

	printf("\nInsert number of printout spots to be displayed for the raw data, usually above 32 for data processing purposes\n");
	fflush(stdin);
	scanf("%d", &SAMPLE_PRINTOUT_SPOTS);

	// Configure WFS camera, use a pre-defined camera resolution
	if((instr.selected_id & DEVICE_OFFSET_WFS10) == 0 && (instr.selected_id & DEVICE_OFFSET_WFS20) == 0) // WFS150/300 instrument
	{   
		printf("\nSelect from available resolutions for WFS150/300:\n");
		printf("\n0 -> 1280x1024 pixels\n");
		printf("\n1 -> 1024x1024 pixels\n");
		printf("\n2 -> 768x768 pixels\n");
		printf("\n3 -> 512x512 pixels\n");
		printf("\n4 -> 320x320 pixels\n");
		fflush(stdin);
		scanf("%d", &device_resolution);
		if (device_resolution == 0)
		{
			printf("\nConfigure WFS camera with resolution index 1280 (1280 x 1024 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_1280, &instr.spots_x, &instr.spots_y))
				handle_errors(err);
		}
		else if (device_resolution == 1)
		{
			printf("\nConfigure WFS camera with resolution index 1024 (1024 x 1024 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_1024, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 2)
		{
			printf("\nConfigure WFS camera with resolution index 768 (768 x 768 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_768, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 3)
		{
			printf("\nConfigure WFS camera with resolution index 512 (512 x 512 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_512, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else
		{
			printf("\nConfigure WFS camera with resolution index 320 (320 x 320 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_320, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}	
	}
	
	if(instr.selected_id & DEVICE_OFFSET_WFS10) // WFS10 instrument
	{
		printf("\nSelect from available resolutions for WFS10:\n");
		printf("\n0 -> 640x480 pixels\n");
		printf("\n1 -> 480x480 pixels\n");
		printf("\n2 -> 360x360 pixels\n");
		printf("\n3 -> 260x260 pixels\n");
		printf("\n4 -> 180x180 pixels\n");
		fflush(stdin);
		scanf("%d", &device_resolution);
		if (device_resolution == 0)
		{
			printf("\nConfigure WFS camera with resolution index 640 (640 x 480 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS10_640, &instr.spots_x, &instr.spots_y))
				handle_errors(err);
		}
		else if (device_resolution == 1)
		{
			printf("\nConfigure WFS camera with resolution index 480 (480 x 480 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS10_480, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 2)
		{
			printf("\nConfigure WFS camera with resolution index 360 (360 x 360 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS10_360, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 3)
		{
			printf("\nConfigure WFS camera with resolution index 260 (260 x 260 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS10_260, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else
		{
			printf("\nConfigure WFS camera with resolution index 180 (180 x 180 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS10_180, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}	
	}
	
	if(instr.selected_id & DEVICE_OFFSET_WFS20) // WFS20 instrument

	{
		printf("\nSelect from available resolutions for WFS20:\n");
		printf("\n0 -> 1440x1080 pixels\n");
		printf("\n1 -> 1080x1080 pixels\n");
		printf("\n2 -> 768x768 pixels\n");
		printf("\n3 -> 512x512 pixels\n");
		printf("\n4 -> 360x360 pixels\n");
		printf("\n5 -> 720x540 pixels, binning 2x2\n");
		printf("\n6 -> 540x540 pixels, binning 2x2\n");
		printf("\n7 -> 384x384 pixels, binning 2x2\n");
		printf("\n8 -> 256x256 pixels, binning 2x2\n");
		printf("\n9 -> 180x180 pixels, binning 2x2\n");
		fflush(stdin);
		scanf("%d", &device_resolution);
		if (device_resolution == 0)
		{
			printf("\nConfigure WFS camera with resolution index 1440 (1440x1080 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_1440, &instr.spots_x, &instr.spots_y))
				handle_errors(err);
		}
		else if (device_resolution == 1)
		{
			printf("\nConfigure WFS camera with resolution index 1080 (1080x1080 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_1080, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 2)
		{
			printf("\nConfigure WFS camera with resolution index 768 (768x768 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_768, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 3)
		{
			printf("\nConfigure WFS camera with resolution index 512 (512x512 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_512, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 4)
		{
			printf("\nConfigure WFS camera with resolution index 360 (360x360 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_360, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}	
		else if (device_resolution == 5)
		{
			printf("\nConfigure WFS camera with resolution index 720_bin (720x540, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_720_BIN2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 6)
		{
			printf("\nConfigure WFS camera with resolution index 540_bin (540x540, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_540_BIN2 , &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 7)
		{
			printf("\nConfigure WFS camera with resolution index 384_bin (384x384, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_384_BIN2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 8)
		{
			printf("\nConfigure WFS camera with resolution index 256_bin (256x256, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_256_BIN2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else
		{
			printf("\nConfigure WFS camera with resolution index 180_bin (180x180, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS20_180_BIN2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}		

	}

	if(instr.selected_id & DEVICE_OFFSET_WFS30) // WFS30 instrument
	{
		printf("\nSelect from available resolutions for WFS30:\n");
		printf("\n0 -> 1936x1216 pixels\n");
		printf("\n1 -> 1216x1216 pixels\n");
		printf("\n2 -> 1024x1024 pixels\n");
		printf("\n3 -> 768x768 pixels\n");
		printf("\n4 -> 512x512 pixels\n");
		printf("\n5 -> 360x360 pixels\n");
		printf("\n6 -> 968x608, binning 2x2\n");
		printf("\n7 -> 608x608, binning 2x2\n");
		printf("\n8 -> 512x512, binning 2x2\n");
		printf("\n9 -> 384x384, binning 2x2\n");
		printf("\n10 -> 256x256, binning 2x2\n");
		printf("\n11 -> 180x180, binning 2x2\n");
		fflush(stdin);
		scanf("%d", &device_resolution);
		if (device_resolution == 0)
		{
			printf("\nConfigure WFS camera with resolution index 1936 (1936x1216 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_1936, &instr.spots_x, &instr.spots_y))
				handle_errors(err);
		}
		else if (device_resolution == 1)
		{
			printf("\nConfigure WFS camera with resolution index 1216 (1216x1216 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_1216, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 2)
		{
			printf("\nConfigure WFS camera with resolution index 1024 (1024x1024 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_1024, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 3)
		{
			printf("\nConfigure WFS camera with resolution index 768 (768x768 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_768, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 4)
		{
			printf("\nConfigure WFS camera with resolution index 512 (512x512 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_512, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}	
		else if (device_resolution == 5)
		{
			printf("\nConfigure WFS camera with resolution index 360 (360x360 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_360, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 6)
		{
			printf("\nConfigure WFS camera with resolution index 968_bin (968x608, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_968_SUB2 , &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 7)
		{
			printf("\nConfigure WFS camera with resolution index 608_bin (608x608, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_608_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 8)
		{
			printf("\nConfigure WFS camera with resolution index 512_bin (512x512, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_512_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 9)
		{
			printf("\nConfigure WFS camera with resolution index 384_bin (384x384, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_384_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 10)
		{
			printf("\nConfigure WFS camera with resolution index 256_bin (256x256, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_256_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else
		{
			printf("\nConfigure WFS camera with resolution index 180_bin (180x180, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS30_180_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}

	}

	if(instr.selected_id & DEVICE_OFFSET_WFS40) // WFS40 instrument
	{
		printf("\nSelect from available resolutions for WFS40:\n");
		printf("\n0 -> 2048x2048 pixels\n");
		printf("\n1 -> 1536x1536 pixels\n");
		printf("\n2 -> 1024x1024 pixels\n");
		printf("\n3 -> 768x768 pixels\n");
		printf("\n4 -> 512x512 pixels\n");
		printf("\n5 -> 360x360 pixels\n");
		printf("\n6 -> 1024x1024, binning 2x2\n");
		printf("\n7 -> 768x768, binning 2x2\n");
		printf("\n8 -> 512x512, binning 2x2\n");
		printf("\n9 -> 384x384, binning 2x2\n");
		printf("\n10 -> 256x256, binning 2x2\n");
		printf("\n11 -> 180x180, binning 2x2\n");
		fflush(stdin);
		scanf("%d", &device_resolution);
		if (device_resolution == 0)
		{
			printf("\nConfigure WFS camera with resolution index 2048 (2048x2048 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_2048, &instr.spots_x, &instr.spots_y))
				handle_errors(err);
		}
		else if (device_resolution == 1)
		{
			printf("\nConfigure WFS camera with resolution index 1536 (1536x1536 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_1536, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 2)
		{
			printf("\nConfigure WFS camera with resolution index 1024 (1024x1024 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_1024, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 3)
		{
			printf("\nConfigure WFS camera with resolution index 768 (768x768 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_768, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 4)
		{
			printf("\nConfigure WFS camera with resolution index 512 (512x512 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_512, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}	
		else if (device_resolution == 5)
		{
			printf("\nConfigure WFS camera with resolution index 360 (360x360 pixels).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_360, &instr.spots_x, &instr.spots_y))
				handle_errors(err);			
		}
		else if (device_resolution == 6)
		{
			printf("\nConfigure WFS camera with resolution index 1024_bin (1024x1024, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_1024_SUB2 , &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 7)
		{
			printf("\nConfigure WFS camera with resolution index 768_bin (768x768, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_768_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 8)
		{
			printf("\nConfigure WFS camera with resolution index 512_bin (512x512, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_512_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 9)
		{
			printf("\nConfigure WFS camera with resolution index 384_bin (384x384, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_384_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else if (device_resolution == 10)
		{
			printf("\nConfigure WFS camera with resolution index 256_bin (256x256, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_256_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
		else
		{
			printf("\nConfigure WFS camera with resolution index 180_bin (180x180, binning 2x2).\n");
			if(err = WFS_ConfigureCam (instr.handle, SAMPLE_PIXEL_FORMAT, CAM_RES_WFS40_180_SUB2, &instr.spots_x, &instr.spots_y))
				handle_errors(err);		
		}
	}
	
	printf("Camera is configured to detect %d x %d lenslet spots.\n\n", instr.spots_x, instr.spots_y);
	
	// set camera exposure time and gain if you don't want to use auto exposure
	// use functions WFS_GetExposureTimeRange, WFS_SetExposureTime, WFS_GetMasterGainRange, WFS_SetMasterGain
	
	// set WFS internal reference plane as the reference to what the centroids will be calculated to
	printf("\nSet WFS to internal reference plane.\n");
	if(err == WFS_SetReferencePlane (instr.handle, SAMPLE_REF_PLANE))
		handle_errors(err);

	if(err == WFS_SetPupil (instr.handle, SAMPLE_PUPIL_CENTROID_X, SAMPLE_PUPIL_CENTROID_Y, SAMPLE_PUPIL_DIAMETER_X, SAMPLE_PUPIL_DIAMETER_Y))
		handle_errors(err);

	
	printf("\nRead camera images:\n");
	
	printf("Image No.     Status     ->   newExposure[ms]   newGainFactor\n");
	
	// do some trials to read a well exposed image
	for(cnt = 0; cnt < SAMPLE_IMAGE_READINGS; cnt++)
	{
		// take a camera image with auto exposure, note that there may several function calls required to get an optimal exposed image
		if(err == WFS_TakeSpotfieldImageAutoExpos (instr.handle, &expos_act, &master_gain_act))
			handle_errors(err);
	
		printf("    %d     ", cnt);
	
		// check instrument status for non-optimal image exposure
		if(err == WFS_GetStatus (instr.handle, &instr.status))
			handle_errors(err);   
	
		if(instr.status & WFS_STATBIT_PTH) printf("Power too high!    ");
		 else
		if(instr.status & WFS_STATBIT_PTL) printf("Power too low!     ");
		 else
		if(instr.status & WFS_STATBIT_HAL) printf("High ambient light!");
		 else
			printf(                                "OK                 ");
		
		printf("     %6.3f          %6.3f\n", expos_act, master_gain_act);
		
		if( !(instr.status & WFS_STATBIT_PTH) && !(instr.status & WFS_STATBIT_PTL) && !(instr.status & WFS_STATBIT_HAL) )
			break; // image well exposed and is usable
	}
	

	// close program if no well exposed image is feasible
	if( (instr.status & WFS_STATBIT_PTH) || (instr.status & WFS_STATBIT_PTL) ||(instr.status & WFS_STATBIT_HAL) )
	{
		printf("\nSample program will be closed because of unusable image quality, press <ENTER>.");
		WFS_close(instr.handle); // required to release allocated driver data
		fflush(stdin);
		getchar();
		exit(1);
	}
	
	
	// get last image (only required to display the image)
	if(err == WFS_GetSpotfieldImage (instr.handle, &ImageBuffer, &rows, &cols))
		handle_errors(err);

	
	// calculate all spot centroid positions using dynamic noise cut option
	if(err == WFS_CalcSpotsCentrDiaIntens (instr.handle, SAMPLE_OPTION_DYN_NOISE_CUT, SAMPLE_OPTION_CALC_SPOT_DIAS))
		handle_errors(err);

	// get centroid result arrays
	if(err == WFS_GetSpotCentroids (instr.handle, *centroid_x, *centroid_y))
		handle_errors(err);


	fp12 = fopen(MLA_FOCAL, "w");
	fprintf(fp12," %5.3f", instr.lenslet_f_um);
    fclose(fp12);
        
	
	// print out some centroid positions
	printf("\nCentroid X Positions in pixels ( %d lenslet elements)\n",SAMPLE_PRINTOUT_SPOTS);
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
			printf(" %8.3f", centroid_x[i][j]);
		printf("\n");  
	}      

	printf("\nCentroid Y Positions in pixels ( %d lenslet elements)\n",SAMPLE_PRINTOUT_SPOTS);
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
			printf(" %8.3f", centroid_y[i][j]);
		printf("\n");  
	}  


	fp6 = fopen(CENTROID_POS_X, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
            {
			fprintf(fp6," %5.3f", centroid_x[i][j]);
			fprintf(fp6,",");
            }  
            fprintf(fp6,"\n"); 
    }
    fclose(fp6);
	
	
	fp7 = fopen(CENTROID_POS_Y, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
            {
			fprintf(fp7," %5.3f", centroid_y[i][j]);
			fprintf(fp7,",");
            }  
            fprintf(fp7,"\n"); 
    }
    fclose(fp7);
    

	printf("\nPress <ENTER> to proceed...");
	getchar();
	

	// get centroid and diameter of the optical beam, you may use this beam data to define a pupil variable in position and size
	// for WFS20: this is based on centroid intensties calculated by WFS_CalcSpotsCentrDiaIntens()
	if(err == WFS_CalcBeamCentroidDia (instr.handle, &beam_centroid_x, &beam_centroid_y, &beam_diameter_x, &beam_diameter_y))
		handle_errors(err);
	
	printf("\nInput beam is measured to:\n");
	printf("Centroid_x = %6.3f mm\n", beam_centroid_x);
	printf("Centroid_y = %6.3f mm\n", beam_centroid_y);
	printf("Diameter_x = %6.3f mm\n", beam_diameter_x);
	printf("Diameter_y = %6.3f mm\n", beam_diameter_y);

	fflush(stdin);
	printf("\nPress <ENTER> to proceed...");
	getchar();
	
	// get pupil centre x,y and diameter, you may use this pupil data to check if the GUI data is in correlation with the saved values
	if(err == WFS_GetPupil (instr.handle, &pupil_centre_x,  &pupil_centre_y, &pupil_diameter_x, &pupil_diameter_y))
		handle_errors(err);
	
	printf("\nPupil is measured to:\n");
	printf("pupil_centre_x = %6.3f mm\n", pupil_centre_x);
	printf("pupil_centre_y = %6.3f mm\n", pupil_centre_y);
	printf("pupil_diameter_x = %6.3f mm\n", pupil_diameter_x);
	printf("pupil_diameter_y = %6.3f mm\n", pupil_diameter_y);

	fflush(stdin);
	printf("\nPress <ENTER> to proceed...");
	getchar();


	// calculate spot deviations to internal reference
	if(err == WFS_CalcSpotToReferenceDeviations (instr.handle, SAMPLE_OPTION_CANCEL_TILT))
		handle_errors(err);
	
	// get spot deviations
	if(WFS_GetSpotDeviations (instr.handle, *deviation_x, *deviation_y))
		handle_errors(err);

    // get spot intensity result arrays
    if(err == WFS_GetSpotIntensities (instr.handle, *intensity))
        handle_errors(err); 
	
	// print out some spot deviations and save to csv file for Matlab Processing
    
    printf("\nCentroid Deviations in X direction have been written to output file %s\n", CENTROID_DEVIATION_FILE_NAME_X);
	printf("\nSpot Deviation X in pixels for %d lenslet elements\n", SAMPLE_PRINTOUT_SPOTS);
    fp1 = fopen(CENTROID_DEVIATION_FILE_NAME_X, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
            {
			fprintf(fp1," %5.3f", deviation_x[i][j]);
			fprintf(fp1,",");
            }  
            fprintf(fp1,"\n"); 
    }
    fclose(fp1);

    printf("\nCentroid Deviations in Y direction have been written to output file %s\n", CENTROID_DEVIATION_FILE_NAME_Y);
	printf("\nSpot Deviation Y in pixels for %d lenslet elements\n", SAMPLE_PRINTOUT_SPOTS);
    fp2 = fopen(CENTROID_DEVIATION_FILE_NAME_Y, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
        {
			fprintf(fp2," %5.3f", deviation_y[i][j]);
		    fprintf(fp2,",");  
        }
        fprintf(fp2,"\n");
	}         
    fclose(fp2);
	
	// Average calculation and printout
	
    fp8 = fopen(CENTROID_DEVIATION_FILE_NAME_X_AVG, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS-1;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS-1;j++)
            {
			deviation_avg_x[i][j] = (deviation_x[i][j] + deviation_x[i][j+1]) / 2.0;
			fprintf(fp8," %5.3f", deviation_avg_x[i][j]);
			fprintf(fp8,",");
            }  
            fprintf(fp8,"\n"); 
    }
    fclose(fp8);

    fp9 = fopen(CENTROID_DEVIATION_FILE_NAME_Y_AVG, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS-1;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS-1;j++)
        {
			deviation_avg_y[i][j] = (deviation_y[i][j] + deviation_y[i][j+1]) / 2.0;
			fprintf(fp9," %5.3f", deviation_avg_y[i][j]);
		    fprintf(fp9,",");  
        }
        fprintf(fp9,"\n");
	}         
    fclose(fp9);

	


    printf("\nField Intensity has been written to output file %s\n", FIELD_INTENSITY_FILE_NAME);
    fp3 = fopen(FIELD_INTENSITY_FILE_NAME, "w");
    for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
    {   
        for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
        {
            fprintf(fp3," %5.3f", intensity[i][j]);
            fprintf(fp3, ",");
        }
        fprintf(fp3,"\n");
    }         
    fclose(fp3);

	fp11 = fopen(FIELD_INTENSITY_FILE_AVG, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS-1;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS-1;j++)
            {
			intensity_avg[i][j] = (intensity[i][j] + intensity[i][j+1]) / 2.0;
			fprintf(fp11," %5.3f", intensity[i][j]);
			fprintf(fp11,",");
            }  
            fprintf(fp11,"\n"); 
    }
    fclose(fp11);

	printf("\nPress <ENTER> to proceed...");
	getchar();
	

	// calculate and printout measured wavefront
	if(err == WFS_CalcWavefront (instr.handle, SAMPLE_WAVEFRONT_TYPE, SAMPLE_OPTION_LIMIT_TO_PUPIL, *wavefront))
		handle_errors(err);
	

	printf("\nWavefront in microns ( %d lenslet elements)\n",SAMPLE_PRINTOUT_SPOTS);
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
			printf(" %8.3f", wavefront[i][j]);
		printf("\n");  
	}     


	fp10 = fopen(WAVEFRONT_ZERNIKE_FIT, "w");
	for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
	{   
		for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
            {
			fprintf(fp10," %5.3f", wavefront[i][j]);
			fprintf(fp10,",");
            }  
            fprintf(fp10,"\n"); 
    }
    fclose(fp10);	
	
	printf("\nPress <ENTER> to proceed...");
	getchar();
	
	
	// calculate wavefront statistics within defined pupil
	if(err == WFS_CalcWavefrontStatistics (instr.handle, &wavefront_min, &wavefront_max, &wavefront_diff, &wavefront_mean, &wavefront_rms, &wavefront_weighted_rms))
		handle_errors(err);
	
	if(err == WFS_GetSpotDiaStatistics (instr.handle, &spot_min, &spot_max, &spot_mean))
		handle_errors(err);
	
	
	printf("\nWavefront Statistics in microns:\n");
	printf("Min          : %8.3f\n", wavefront_min);
	printf("Max          : %8.3f\n", wavefront_max);
	printf("Diff         : %8.3f\n", wavefront_diff);
	printf("Mean         : %8.3f\n", wavefront_mean);
	printf("RMS          : %8.3f\n", wavefront_rms);
	printf("Weigthed RMS : %8.3f\n", wavefront_weighted_rms);
	printf("Minimum spot diameter:  %8.3f\n", spot_min);
	printf("Maximum spot diameter:  %8.3f\n", spot_max);
	printf("Mean average of spot diameters:  %8.3f\n", spot_mean);
	
	
	printf("\nPress <ENTER> to proceed...");
	getchar();
	
	if(err == WFS_GetSpotReferencePositions (instr.handle, &spot_ref_x, &spot_ref_y))
		handle_errors(err);
	
	
	printf("\nX reference possitions in pixels has been written to output file %s\n",REFERENCE_SPOT_POS_X );
    fp4 = fopen(REFERENCE_SPOT_POS_X , "w");
    for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
    {   
        for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
        {
            fprintf(fp4," %5.3f", spot_ref_x[i][j]);
            fprintf(fp4, ",");
        }
        fprintf(fp4,"\n");
    }         
    fclose(fp4);
	
	
	printf("\nY reference possitions in pixels has been written to output file %s\n", REFERENCE_SPOT_POS_Y );
    fp5 = fopen(REFERENCE_SPOT_POS_Y, "w");
    for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
    {   
        for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
        {
            fprintf(fp5," %5.3f", spot_ref_y[i][j]);
            fprintf(fp5, ",");
        }
        fprintf(fp5,"\n");
    }         
    fclose(fp5);
	
	
	printf("\nPress <ENTER> to proceed...");
	getchar();
	
	// calculate Zernike coefficients
	printf("\nZernike fit up to order %d:\n",SAMPLE_ZERNIKE_ORDERS);
	zernike_order = SAMPLE_ZERNIKE_ORDERS; // pass 0 to function for auto Zernike order, choosen order is returned
	if(err == WFS_ZernikeLsf (instr.handle, &zernike_order, zernike_um, zernike_orders_rms_um, &roc_mm)) // calculates also deviation from centroid data for wavefront integration
		handle_errors(err);
		
	printf("\nZernike Mode    Coefficient\n");
	for(i=0; i < zernike_modes[SAMPLE_ZERNIKE_ORDERS]; i++)
	{
		printf("  %2d         %9.3f\n",i, zernike_um[i]);
	}

	
	
	printf("\nEnter measurement loop with output to file 0/1?");
	fflush(stdin);
	selection = getchar() - '0';

	if(selection == 1)
	{
		printf("\nMeasurement data is continuously written into file %s.\n", OUTPUT_FILE_NAME);
	
		printf("\nInsert number of iterations the loop will run for\n");
		fflush(stdin);
		scanf("%d", &iterations);
		for (int i = 0 ; i<iterations; i++)
		{   
			// take a camera image with auto exposure
			if(err == WFS_TakeSpotfieldImageAutoExpos (instr.handle, &expos_act, &master_gain_act))
				handle_errors(err);
	
			// calculate all spot centroid positions using dynamic noise cut option
			if(err == WFS_CalcSpotsCentrDiaIntens (instr.handle, SAMPLE_OPTION_DYN_NOISE_CUT, SAMPLE_OPTION_CALC_SPOT_DIAS))
				handle_errors(err);

			// calculate spot deviations to internal reference
			if(err == WFS_CalcSpotToReferenceDeviations (instr.handle, SAMPLE_OPTION_CANCEL_TILT))
				handle_errors(err);
	

			// calculate measured wavefront
			if(err == WFS_CalcWavefront (instr.handle, SAMPLE_WAVEFRONT_TYPE, SAMPLE_OPTION_LIMIT_TO_PUPIL, *wavefront))
				handle_errors(err);
	
			// calculate wavefront statistics within defined pupil
			if(err == WFS_CalcWavefrontStatistics (instr.handle, &wavefront_min, &wavefront_max, &wavefront_diff, &wavefront_mean, &wavefront_rms, &wavefront_weighted_rms))
				handle_errors(err);
	
			// calculate Zernike coefficients
			zernike_order = SAMPLE_ZERNIKE_ORDERS; // pass 0 to function for auto Zernike order, choosen order is returned
			if(err == WFS_ZernikeLsf (instr.handle, &zernike_order, zernike_um, zernike_orders_rms_um, &roc_mm)) // calculates also deviation from centroid data for wavefront integration
				handle_errors(err);
		
	
			// copy some values into a text file, overwrite old file content
			if((fp = fopen (OUTPUT_FILE_NAME, "w")) != NULL)
			{   
				fprintf(fp, "Wavefront results in um:\n");
				fprintf(fp, "%s %8.3f\n", "PV   ", wavefront_diff);
				fprintf(fp, "%s %8.3f\n", "RMS  ", wavefront_rms);
		
				fprintf(fp, "\nZernike amplitudes in um:\n");
				for(i=0;i < zernike_modes[SAMPLE_ZERNIKE_ORDERS]; i++)
					fprintf(fp, "%2d    %8.3f\n", i, zernike_um[i]);

				fclose(fp);
			}
			
		}

	}
	
	
	// enter highspeed mode for WFS10 and WFS20 instruments only
	if((instr.selected_id & DEVICE_OFFSET_WFS10) || (instr.selected_id & DEVICE_OFFSET_WFS20)) // WFS10 or WFS20 instrument
	{
		printf("\nEnter Highspeed Mode 0/1?");
		fflush(stdin);
		selection = getchar() - '0';
	
		if(selection == 1)
		{   
			if(err == WFS_SetHighspeedMode (instr.handle, SAMPLE_OPTION_HIGHSPEED, SAMPLE_OPTION_HS_ADAPT_CENTR, SAMPLE_HS_NOISE_LEVEL, SAMPLE_HS_ALLOW_AUTOEXPOS))
				handle_errors(err);
						  
			if(err == WFS_GetHighspeedWindows (instr.handle, &hs_win_count_x, &hs_win_count_y, &hs_win_size_x, &hs_win_size_y, hs_win_start_x, hs_win_start_y))handle_errors(err);
	
			printf("\nCentroid detection windows are defined as follows:\n"); // refere to WFS_GetHighspeedWindows function help
			printf("Count_x = %3d, Count_y = %3d\n", hs_win_count_x, hs_win_count_y);
			printf("Size_x  = %3d, Size_y  = %3d\n", hs_win_size_x, hs_win_size_y);
			printf("Start coordinates x: ");
			for(i=0;i<hs_win_count_x;i++)
				printf("%3d ", hs_win_start_x[i]);
			printf("\n");
			printf("Start coordinates y: ");
			for(i=0;i<hs_win_count_y;i++)
				printf("%3d ", hs_win_start_y[i]);
			printf("\n");
			
			
			fflush(stdin);
			printf("\nPress <ENTER> to proceed...");
			getchar();
			
		
		    // take a camera image with auto exposure, this is also supported in highspeed-mode
		    if(err == WFS_TakeSpotfieldImageAutoExpos (instr.handle, &expos_act, &master_gain_act))
			    handle_errors(err);

		    printf("\nexposure = %6.3f ms, gain =  %6.3f\n", expos_act, master_gain_act);
			

			// get centroid and diameter of the optical beam, these data are based on the detected centroids
			if(err == WFS_CalcBeamCentroidDia (instr.handle, &beam_centroid_x, &beam_centroid_y, &beam_diameter_x, &beam_diameter_y))
				handle_errors(err);
	
			printf("\nInput beam is measured to:\n");
			printf("Centroid_x = %6.3f mm\n", beam_centroid_x);
			printf("Centroid_y = %6.3f mm\n", beam_centroid_y);
			printf("Diameter_x = %6.3f mm\n", beam_diameter_x);
			printf("Diameter_y = %6.3f mm\n", beam_diameter_y);

			fflush(stdin);
			printf("\nPress <ENTER> to proceed...");
			getchar();
			
			// Info: calling WFS_CalcSpotsCentrDiaIntens() is not required because the WFS10/WFS20 camera itself already did the calculation
			
			// get centroid result arrays
			if(err == WFS_GetSpotCentroids (instr.handle, *centroid_x, *centroid_y))
				handle_errors(err);
			
			// print out some centroid positions
			printf("\nCentroid X Positions in pixels ( %d lenslet elements)\n",SAMPLE_PRINTOUT_SPOTS);
			for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
			{   
				for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
					printf(" %8.3f", centroid_x[i][j]);
				printf("\n");  
			}      

			printf("\nCentroid Y Positions in pixels ( %d lenslet elements)\n",SAMPLE_PRINTOUT_SPOTS);
			for(i=0;i<SAMPLE_PRINTOUT_SPOTS;i++)
			{   
				for(j=0;j<SAMPLE_PRINTOUT_SPOTS;j++)
					printf(" %8.3f", centroid_y[i][j]);
				printf("\n");  
			}      		
			printf("\nThe following wavefront and Zernike calculations can be done identical to normal mode.\n");
		}
	}

	ROWS = SAMPLE_PRINTOUT_SPOTS;
	COLS = SAMPLE_PRINTOUT_SPOTS;
	printf("\n Insert number of columns and rows to be deleted from data, due to NaN or QNaN values \n");
	fflush(stdin);
	scanf("%d", &rows_columns_del);
	printf("\n This deletes the first and the last rows/columns %d amount of times\n",rows_columns_del);
	printf("\n Press <ENTER> to proceed...\n");
	getchar();

	if (rows_columns_del == 0)
	{
		printf("\n No rows and columns have been deleted, end of program \n");
		WFS_close(instr.handle);
	}
	else
	{
		while(1)
		{
			for (int i = 0; i <rows_columns_del ; i++)
			{
				// Deleting first row and first column
			    deleteRow(centroid_x, 1, ROWS, COLS);	
			    deleteRow(centroid_y, 1, ROWS, COLS);
			    deleteRow(spot_ref_y, 1, ROWS, COLS);
			    deleteRow(spot_ref_x, 1, ROWS, COLS);
			    deleteRow(intensity, 1, ROWS, COLS);
				deleteRow(deviation_y, 1, ROWS, COLS);
			    deleteRow(deviation_x, 1, ROWS, COLS);
				deleteRow(wavefront, 1, ROWS, COLS);
				deleteRow(intensity_avg, 1 , ROWS, COLS);
				deleteRow(deviation_avg_x, 1, ROWS, COLS);
				deleteRow(deviation_avg_y, 1, ROWS, COLS);
			    ROWS--;

			    deleteCol(centroid_x, 1, ROWS, COLS);	
			    deleteCol(centroid_y, 1, ROWS, COLS);
			    deleteCol(spot_ref_y, 1, ROWS, COLS);
			    deleteCol(spot_ref_x, 1, ROWS, COLS);
			    deleteCol(intensity, 1, ROWS, COLS);
				deleteCol(deviation_y, 1, ROWS, COLS);
			    deleteCol(deviation_x, 1, ROWS, COLS);
				deleteCol(wavefront, 1, ROWS, COLS);
				deleteCol(intensity_avg, 1, ROWS, COLS);
				deleteCol(deviation_avg_x, 1, ROWS, COLS);
				deleteCol(deviation_avg_y, 1, ROWS, COLS);
			    COLS--;

			    //Deleting last row and last column
			    deleteRow(centroid_x, ROWS, ROWS, COLS);	
			    deleteRow(centroid_y, ROWS, ROWS, COLS);
			    deleteRow(spot_ref_y, ROWS, ROWS, COLS);
			    deleteRow(spot_ref_x, ROWS, ROWS, COLS);
			    deleteRow(intensity, ROWS, ROWS, COLS);
				deleteRow(deviation_y, ROWS, ROWS, COLS);
			    deleteRow(deviation_x, ROWS, ROWS, COLS);
				deleteRow(wavefront, ROWS, ROWS, COLS);
				deleteRow(intensity_avg, ROWS, ROWS, COLS);
				deleteRow(deviation_avg_x, ROWS, ROWS, COLS);
				deleteRow(deviation_avg_y, ROWS, ROWS, COLS);
			    ROWS--;

			    deleteCol(centroid_x, COLS, ROWS, COLS);	
			    deleteCol(centroid_y, COLS, ROWS, COLS);
			    deleteCol(spot_ref_y, COLS, ROWS, COLS);
			    deleteCol(spot_ref_x, COLS, ROWS, COLS);
			    deleteCol(intensity, COLS, ROWS, COLS);
				deleteCol(deviation_y, COLS, ROWS, COLS);
			    deleteCol(deviation_x, COLS, ROWS, COLS);
				deleteCol(wavefront, COLS, ROWS, COLS);
				deleteCol(intensity_avg, COLS, ROWS, COLS);
				deleteCol(deviation_avg_x, COLS, ROWS, COLS);
				deleteCol(deviation_avg_y, COLS, ROWS, COLS);
			    COLS--;
			}

			// Writing truncated data to files

			fp6 = fopen(PROCESS_CENTROID_POS_X, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		            {
					fprintf(fp6," %5.3f", centroid_x[i][j]);
					fprintf(fp6,",");
		            }  
		            fprintf(fp6,"\n"); 
		    }
		    fclose(fp6);
			
			
			fp7 = fopen(PROCESS_CENTROID_POS_Y, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		            {
					fprintf(fp7," %5.3f", centroid_y[i][j]);
					fprintf(fp7,",");
		            }  
		            fprintf(fp7,"\n"); 
		    }
		    fclose(fp7);

		    fp4 = fopen(PROCESS_REFERENCE_SPOT_POS_X , "w");
		    for(i=0;i<COLS;i++)
		    {   
		        for(j=0;j<COLS;j++)
		        {
		            fprintf(fp4," %5.3f", spot_ref_x[i][j]);
		            fprintf(fp4, ",");
		        }
		        fprintf(fp4,"\n");
		    }         
		    fclose(fp4);

	    	fp11 = fopen(PROCESS_FIELD_INTENSITY_FILE_AVG, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		            {
					fprintf(fp11," %5.3f", intensity_avg[i][j]);
					fprintf(fp11,",");
		            }  
		            fprintf(fp11,"\n"); 
		    }
		    fclose(fp11);
			
			
		    fp5 = fopen(PROCESS_REFERENCE_SPOT_POS_Y, "w");
		    for(i=0;i<COLS;i++)
		    {   
		        for(j=0;j<COLS;j++)
		        {
		            fprintf(fp5," %5.3f", spot_ref_y[i][j]);
		            fprintf(fp5, ",");
		        }
		        fprintf(fp5,"\n");
		    }         
		    fclose(fp5);

		    fp1 = fopen(PROCESS_CENTROID_DEVIATION_FILE_NAME_X, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		            {
					fprintf(fp1," %5.3f", deviation_x[i][j]);
					fprintf(fp1,",");
		            }  
		            fprintf(fp1,"\n"); 
		    }
		    fclose(fp1);

		    fp2 = fopen(PROCESS_CENTROID_DEVIATION_FILE_NAME_Y, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		        {
					fprintf(fp2," %5.3f", deviation_y[i][j]);
				    fprintf(fp2,",");  
		        }
		        fprintf(fp2,"\n");
			}         
		    fclose(fp2); 

		    fp3 = fopen(PROCESS_FIELD_INTENSITY_FILE_NAME, "w");
		    for(i=0;i<COLS;i++)
		    {   
		        for(j=0;j<COLS;j++)
		        {
		            fprintf(fp3," %5.3f", intensity[i][j]);
		            fprintf(fp3, ",");
		        }
		        fprintf(fp3,"\n");
		    }         
		    fclose(fp3);

			fp8 = fopen(PROCESS_CENTROID_DEVIATION_FILE_NAME_X_AVG, "w");
			for(i=0;i<COLS-1;i++)
			{   
				for(j=0;j<COLS-1;j++)
				{
					fprintf(fp8," %5.3f", deviation_avg_x[i][j]);
					fprintf(fp8,",");
				}  
				fprintf(fp8,"\n"); 
			}
			fclose(fp8);

			fp9 = fopen(PROCESS_CENTROID_DEVIATION_FILE_NAME_Y_AVG, "w");
			for(i=0;i<COLS-1;i++)
			{   
				for(j=0;j<COLS-1;j++)
				{
					fprintf(fp9," %5.3f", deviation_avg_y[i][j]);
					fprintf(fp9,",");  
				}
				fprintf(fp9,"\n");
			}         
			fclose(fp9);

			fp10 = fopen(PROCESS_WAVEFRONT_ZERNIKE_FIT, "w");
			for(i=0;i<COLS;i++)
			{   
				for(j=0;j<COLS;j++)
		            {
					fprintf(fp10," %5.3f", wavefront[i][j]);
					fprintf(fp10,",");
		            }  
		            fprintf(fp10,"\n"); 
		    }
		    fclose(fp10);

	    	printf("\n New truncated data by %d rows and columns has been written to PROCESS_ files \n", rows_columns_del);
	    	printf("\n If you wish continue truncating rows and columns, based on visual inspection of the .CSV files, please insert a new number, else 0 if you wish to stop\n");
			fflush(stdin);
			scanf("%d", &rows_columns_del);
			printf("\n This deletes the first and the last rows/columns %d amount of times\n",rows_columns_del);
			printf("\n Press <ENTER> to proceed...\n");
			getchar();
			if (rows_columns_del == 0)
				break;
			else
				continue;
		}
	printf("\nEnd of program, press <ENTER> to exit.");
	getchar();
	WFS_close(instr.handle);
	}
}


// Delete row function

void deleteRow(float array[], int row, int rows, int cols)
{
    int newsize, index, i;
    newsize = (rows - 1) * cols;
    if (row != rows){
        index = cols * row;
        for (i = index; i<newsize; i++){
            array[i] = array[i+cols];
        }
    }
}

// Delete column function

void deleteCol(float array[], int col, int rows, int cols)
{
    int newsize, index, i, j;
    newsize = rows * (cols - 1);
    index = col;
    for (i = index; i<newsize; i++){
        int di=1 +(1-index)/(cols-1);
        array[i] = array[i+di];
    }
}


/*===============================================================================================================================
  Handle Errors
  This function retrieves the appropriate text to the given error number and closes the connection in case of an error
===============================================================================================================================*/
void handle_errors (int err)
{
	char buf[WFS_ERR_DESCR_BUFFER_SIZE];

	if(!err) return;

	// Get error string
	WFS_error_message (instr.handle, err, buf);

	if(err < 0) // errors
	{
		printf("\nWavefront Sensor Error: %s\n", buf);

		// close instrument after an error has occured
		printf("\nSample program will be closed because of the occured error, press <ENTER>.");
		WFS_close(instr.handle); // required to release allocated driver data
		fflush(stdin);
		getchar();
		exit(1);
	}
}



/*===============================================================================================================================
	Select Instrument
===============================================================================================================================*/
int select_instrument (int *selection, ViChar resourceName[])
{
	int            i,err=0,instr_cnt;
	ViInt32        device_id;
	int            in_use;
	char           instr_name[WFS_BUFFER_SIZE];
	char           serNr[WFS_BUFFER_SIZE];
	char           strg[WFS_BUFFER_SIZE];

	// Find available instruments
	if(err == WFS_GetInstrumentListLen (VI_NULL, &instr_cnt))
		handle_errors(err);
		
	if(instr_cnt == 0)
	{
		printf("No Wavefront Sensor instrument found!\n");
		return 0;
	}

	// List available instruments
	printf("Available Wavefront Sensor instruments:\n\n");
	
	for(i=0;i<instr_cnt;i++)
	{
		if(err == WFS_GetInstrumentListInfo (VI_NULL, i, &device_id, &in_use, instr_name, serNr, resourceName))
			handle_errors(err);
		
		printf("%3d   %s    %s    %s\n", device_id, instr_name, serNr, (!in_use) ? "" : "(inUse)");
	}

	// Select instrument
	printf("\nSelect a Wavefront Sensor instrument: ");
	fflush(stdin);
	
	fgets (strg, WFS_BUFFER_SIZE, stdin);
	*selection = atoi(strg);

	// get selected resource name
	for(i=0;i<instr_cnt;i++)
	{   
		if(err == WFS_GetInstrumentListInfo (VI_NULL, i, &device_id, &in_use, instr_name, serNr, resourceName))
		   handle_errors(err);
		
		if(device_id == *selection)
			break; // resourceName fits to device_id
	}
	
	return *selection;
}


/*===============================================================================================================================
	Select MLA
===============================================================================================================================*/
int select_mla (int *selection)
{
	int            i,err=0,mla_cnt;

	// Read out number of available Microlens Arrays 
	if(err == WFS_GetMlaCount (instr.handle, &instr.mla_cnt))
		handle_errors(err);

	// List available Microlens Arrays
	printf("\nAvailable Microlens Arrays:\n\n");
	for(i=0;i<instr.mla_cnt;i++)
	{   
		if(WFS_GetMlaData (instr.handle, i, instr.mla_name, &instr.cam_pitch_um, &instr.lenslet_pitch_um, &instr.center_spot_offset_x, &instr.center_spot_offset_y, &instr.lenslet_f_um, &instr.grd_corr_0, &instr.grd_corr_45))
			handle_errors(err);   
	
		printf("%2d  %s   CamPitch=%6.3f LensletPitch=%8.3f\n", i, instr.mla_name, instr.cam_pitch_um, instr.lenslet_pitch_um);
	}
	
	// Select MLA
	printf("\nSelect a Microlens Array: ");
	fflush(stdin);
	*selection = getchar() - '0';
	if(*selection < -1)
		*selection = -1; // nothing selected

	return *selection;
}


/*===============================================================================================================================
	End of source file
===============================================================================================================================*/
