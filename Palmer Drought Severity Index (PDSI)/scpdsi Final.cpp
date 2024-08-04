#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <ctype.h>
#include <direct.h>  
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//			UNIVERSIDAD DISTRITAL FRANCISCO JOSÉ DE CALDAS-BOGOTÁ-COLOMBIA
// 							AUTOR : CARLOS MÉNDEZ

// Contacto : camendezv@correo.udistrital.edu.co
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

// Créditos, versiones anteriores de referencia 
// University of Nebraska - Lincoln            Jul 15 2003 
// Rob Reed and Nate Wells
// Tom Heddinghaus en NCEP
// programas anteriores 
// scpdsi.cpp
// weekly_pdsi.cpp
// ----------------------------------------------------------------------------------------
// 					CARACTERÍSTICAS DEL PROGRAMA 
//************** Calcula el Índice de severidad de sequía de Palmer (PDSI) 
//************** Calcula el Índice de sequía hidrológica de Palmer (PHDI)
//************** Calcula el Índice de severidad de sequía de Palmer modificado o ponderado(WPLM)
//************** Calcula el Índice Z de Palmer (ZIND)
//************** El valor de los índices se calculan en escalas de tiempo mensuales
//************** Acepta un formato de texto con los parámetros necesarios:
//************** Utiliza el método de Evapotranspiración de Thornthwaite 
// ----------------------------------------------------------------------------------------
//************** PARÁMETROS DE ENTRADA

//************** VARIABLE DE PRECIPITACIÓN  (unidades en milímetros) 
//************** Precipitación Total Mensual
//************** Nombre del archivo ------>  monthly_P

//************** VARIABLE TEMETPRATURA (unidades en grados centigrados)
//************** Temperatura Media Mensual 
//************** Nombre del archivo ------>  monthly_T

//************** VARIABLE CAPACIDAD DE ALMACENAMIENTO/RETENCIÓN/APROVECHAMIENTO DE AGUA O DE HÚMEDAD DEL SUELO (unidades de milímetros) 
//************** AWC (capacidad de retencion de agua disponible ) y Latitud Norte o Sur:
//************** Nombre del archivo ------> parameter

//************** VARIABLE EVAPOTRANSPIRACIÓN POTENCIAL(unidades de milímetros)
//************** EVAPOTRANSPIRACIÓN POTENCIAL ajustada o corregida con los valores de latitud, Horas de sol y Temperatura
//************** Nombre del archivo ------> monthly_ETP


// Se modificó el método de calibración para calibrar la característica climática según
// en los cuartiles en lugar del máximo y mínimo
// ------------------------------------------------ -----------------------------
// ------ Archivos de salida:
//
// Los archivos de salida de la tabla terminan con .tbl 
// los archivos de salida de la columna terminan en .clm. 
//
// PDSI.tbl y / o PDSI.clm:
//   Los valores del índice de severidad de la sequía de Palmer
//
// PHDI.tbl y / o PHDI.clm:
//   Valores del índice de sequía hidrológica de Palmer.
//
// WPLM.tbl y / o WPLM.clm:
//   El índice de Palmer "ponderado". Un promedio de X1 y X3, o X1 y X2
//   ponderado por la probabilidad de que finalice el hechizo actual. 
// 
// ZIND.tbl y / o ZIND.clm:
//   El índice Z, o la anomalía de la humedad
//
// WB.tbl
//   Los coeficientes de balance de agua (Alfa, Beta, Gamma y Delta) para cada
//   semana / mes.
//
// bigTable.tbl
//   Z,% Prob (final), X1, X2 y X3 para cada mes.
//
// potenciales
//   P, ETP, PR, PRO, PL y P - ETP para cada mes.
//
// dvalue
// La salida de humedad, d, para cada mes.
//
// Esto define el número de tipo como doble. Esto se usa para cambiar fácilmente
// los tipos de variables del PDSI.

// Esto define el número de tipo como doble. Esto se usa para cambiar fácilmente
// los tipos de variables del PDSI.

typedef double number;
typedef int flag;
#define min(a,b) ((a) < (b) ? (a) : (b));
#define MISSING -99.00
struct node {            // A structure for a node
public:
  number key;            // Where the data is stored
  struct node *next;     // Where the next node is
  struct node *previous; // Where the previous node is
};
class llist {            // A linked list class

private:
  node *head;            // A pointer to the head of the linked-list
  int size;
  number kthLargest(int k);  //returns the kth largest number in the list
  number long_select(int k); //returns the kth largest number in the list
  number percentile(double percentage); //returns the specified percentile
  
public:
  llist();               // The constructor
  ~llist();              // The destructor
  void insert(number x); 
  int is_empty();// Returns 1 if the llist is empty 0 otherwise
  int get_size();
  void sumlist(number &prev_sum, int sign);//sums items in list into prev_sum
  number head_remove();  // remove the first node and returns its key
  number tail_remove();  // remove the last node and returns its key
  number sumlist();  // Sums the items in list
  number maxlist();
  number minlist();
  number quartile(int q);        //returns the qth quartile
  number safe_percentile(double percentage); //safe version 
  number* returnArray();
  node *set_node(node *set=NULL, number x=0);
  friend void copy(llist &L1,const llist &L2); // Copies L2 to L1
};
class pdsi {

public:
  
  pdsi();
  ~pdsi();
  
  void set_flags(int num_flags, char *flags[]);
  void MonthlyPDSI();
  void SCMonthlyPDSI();
  void Write();
  void Write(char* dir);
  
  number getPDSI(int period, int year);
  number getZIND(int period, int year);
  number getPHDI(int period, int year);
  number getWPLM(int period, int year);
  
  number* getPDSIArray(int &size);
  number* getPDSIArray(int start_per, int start_yr, int end_per, int end_yr, int &size);
  number* getZINDArray(int &size);
  number* getZINDArray(int start_per, int start_yr,int end_per, int end_yr, int &size);
  number* getWPLMArray(int &size);
  number* getWPLMArray(int start_per, int start_yr, int end_per, int end_yr, int &size);
  number* getPHDIArray(int &size);
  number* getPHDIArray(int start_per, int start_yr,int end_per, int end_yr, int &size);
  number* getYearArray(int &size);
  number* getPerArray(int &size);     
       
private:
	
  bool Weekly;
  bool Monthly;
  bool SCMonthly;
  
  int startyear;
  int endyear;
  int totalyears;
  int period_length;        //set to 1 for monthly, otherwise, legth of period
  int num_of_periods;       //number of periods of period_length in a year.
  
  flag bug;
  flag output_mode;
  flag verbose;
  flag s_year;
  flag e_year;
  flag extra;
  flag metric;
  flag south;
  flag nadss;
  flag setCalibrationStartYear;
  flag setCalibrationEndYear;
  
  int calibrationStartYear;
  int calibrationEndYear;
  int currentCalibrationStartYear;
  int currentCalibrationEndYear;
  int nStartYearsToSkip;
  int nEndYearsToSkip;
  int nCalibrationYears;
  int nStartPeriodsToSkip;
  int nEndPeriodsToSkip;
  int nCalibrationPeriods;
  
  number TLA; // The negative tangent of latitude is used in calculating PE
  number AWC; // The soils water capacity
  number tolerance; // The tolerance for various comparisons
  number TNorm[52];
  
  number T[52];
  number P[52];
  number ET;            // Actual evapotranspiration
  number R;             // Actual soil recharge 
  number L;             // Actual loss
  number RO;            // Actual runoff
  number PE;            // Potential evapotranspiration
  number PR;            // Potential soil recharge
  number PL;            // Potential Loss
  number PRO;           // Potential runoff
  number Su;            // Underlying soil moisture
  number Ss;            // Surface soil moisture
  
  // Arrays to calculate the values of the potentials
  number ETSum[52];
  number RSum[52];
  number LSum[52];
  number ROSum[52];
  number PESum[52];
  number PRSum[52];
  number PLSum[52];
  number PROSum[52];
  number PSum[52];
  
 // constants for convert the potentials values to real values
  number Alpha[52];
  number Beta[52];
  number Gamma[52];
  number Delta[52];
  
  
  number Phat;
  
  number d;     // Departure from normal for a period
  number D[52]; // Sum of the absolute value of all d values by period
  number k[52]; // Palmer's k' constant by period
  number K;     // The final K value for a period
  number Z;     // The z index for a period (Z=K*d)
  
// constants for adjust the PDSI
// var m ----->   slope of linear regression
// var b ----->   intercept of linear regression
  number drym;  // dry duration m
  number dryb;  // dry duration b 
  number wetm;  // wet duration m
  number wetb;  // wet duration b
  number dry_ratio; // Dry ratio
  number wet_ratio; // Wet/humidity ratio 
  
  number X1;    // Wet index for a month/week
  number X2;    // Dry index for a month/week
  number X3;    // Index for an established wet or dry spell
  

  number X;     // Current period's pdsi value before backtracking  
  number Prob;  // Prob=V/Q*100
  number V;     // Sumation of effective wetness/dryness
  number Q;     // Z needed for an end plus last period's V
 
  number DSSqr[52];
  number DEPSum[52];
  number DKSum;
  number SD;
  number SD2;
  
  llist Xlist;//final list of PDSI values
  llist altX1;//list of X1 values
  llist altX2;//list of X2 values
  llist XL1;
  llist XL2;
  llist XL3;
  llist ProbL;
  llist ZIND;
  llist PeriodList;
  llist YearList;
  
  char input_dir[128];
  char output_dir[128];
  
  int initialize();
  int check_input(FILE *in);
  int GetTemp(FILE * In, number *A, int max); // Gets the Temp
  int GetPrecip(FILE *In, number *A, int max); // Gets the Precip
  int ReprintFile(char* src, char *des);
  
  void GetParam(FILE * Param);          // Gets Paramaters Su and TLA
  void CalcMonPE(int month, int year);
  void CalcPR();// Calculates Potential Recharge
  void CalcPL();// Calculates Potential Loss
  void CalcPRO();// Calculates Potential Runoff
  void CalcActual(int per);// Calculates Actual values
  void SumAll(); // Creates sums of Actual and Potential values 
  void Calcd();     // Finds the values for d, phat and then finds D
  void CalcK();     // Calculates K
  void CalcOrigK();  // This function does it all. It computes K, the Z-index 
  void CalcZ();     // Calculates the Z-index
  void CalcX();     // Calculates the PDSI and X1, X2, and X3
  void CalcOneX(FILE* table, int period_number, int year);
  void Calibrate(); // Calibrates the PDSI
  void Backtrack(number X1, number X2);
  void ChooseX(number& newX, number& newX1, number& newX2, number& newX3, int bug);
  
  void CalcWBCoef(); // this function print the constants Alpha, Beta, Gamma and Delta
  void CalcDurFact(number &slope, number &intercept, int sign);
  void LeastSquares(int *x, number *y, int n, int sign, number &slope, number &intercept); 
  void FinalWrite(char* dir);
  void MoveFiles(char *dir);
  
  number get_Z_sum(int length, int sign);
  number CalcThornA(number I);// Calculates Thornthwaite Exponent
  number getValue(llist &List, int period, int year);
  number* getSubArray(llist &List, int start_per, int start_yr, int end_per, int end_yr, int &size);
  
  inline int is_int(char *string,int length);
  inline int is_flt(char *string,int length);
  
};

int dir_exists(char *dir);
int create_dir(char *path);
void select(number a[], int l, int r, int k);
int partition(number a[], int left, int right);
void exch(number &x, number &y);
int numEntries(FILE *in);

int main(int argc,char *argv[]) {
  pdsi PDSI;
  PDSI.set_flags(argc,argv); // Sets the flags of PDSI 
  PDSI.MonthlyPDSI();
  PDSI.Write("monthly/original");
  PDSI.SCMonthlyPDSI();
  PDSI.Write("monthly/self_cal");
  return 0;
}
pdsi::pdsi() {
  strcpy(input_dir,"./");
  strcpy(output_dir,"./");
  period_length = 1;
  num_of_periods = 52;
  verbose=1;
  bug=0;
  output_mode=0;
  tolerance=0.00001;
  metric=0;
  nadss=0;
  setCalibrationStartYear=0;
  setCalibrationEndYear=0;
}
pdsi::~pdsi() {
  if(verbose > 0)
    printf("\   * Calculations Complete\n");
}
int pdsi::check_input(FILE *in) {
  float temp;
  int count = 0;
  int min_years = 25;
  if(in == NULL)
    return -1;
  while(fscanf(in, "%f",&temp) != EOF)
  {	
    //if temp is either MISSING, or a year, don't count it.
    if(temp != MISSING && temp <= 999)
      count++;
  }  
    if(count < (min_years * 12) )
      return 0;
    else 
      return 1;
}
int pdsi::initialize() {
  char filename[170], base[170];
  FILE* test;
  while(!Xlist.is_empty())
    Xlist.tail_remove();
  while(!altX1.is_empty())
    altX1.tail_remove();
  while(!altX2.is_empty())
    altX2.tail_remove();
  while(!XL1.is_empty())
    XL1.tail_remove();
  while(!XL2.is_empty())
    XL2.tail_remove();
  while(!XL3.is_empty())
    XL3.tail_remove();
  while(!ProbL.is_empty())
    ProbL.tail_remove();
  while(!ZIND.is_empty())
    ZIND.tail_remove();
  while(!PeriodList.is_empty())
    PeriodList.tail_remove();
  while(!YearList.is_empty())
    YearList.tail_remove();
  //check to make sure the necessary input files exist
  if(strlen(input_dir)>1)
    strcpy(base,input_dir);
  else
    strcpy(base,"./");
  if(Monthly || SCMonthly)
  {
    strcpy(filename,base);
    strcat(filename,"parameter");
    strcpy(filename,base);
    strcat(filename,"monthly_P");
    strcpy(filename,base);
    strcat(filename,"monthly_T");
  }
  return 1;
}
void pdsi::set_flags(int num_flags,char *flags[]) {
  int read_from =1;
  unsigned int n;
  flag year=0, both=0, out=0, e_flag=0, t_year=0;
  flag in_dir = -1, out_dir = -1;
  s_year=0;
  e_year=0;
  float t;
  FILE * scn;
  char filename[170];
  if(strlen(input_dir)>1)
  { 	
    strcpy(filename,input_dir);
  }
  //**********************************
  if((scn=fopen(filename,"r"))==NULL) 
  {
	    if(strlen(input_dir)>1)
			{
	      strcpy(filename,input_dir);
	      strcat(filename,"monthly_T");
	    	}
	    else
		  strcpy(filename,"monthly_T");
		  strcpy(filename,input_dir);
		  strcat(filename,"monthly_T");
	      	if((scn=fopen(filename,"r"))==NULL) 
			{    
	     	exit(1);
    		}
	      //scan monthly_T for start and end year
	      fscanf(scn,"%d",&startyear);
		      while((fscanf(scn,"%f",&t))!=EOF)
			  {
				if(t>999)
				  endyear=(int)t;
	      	  }
  }
  //**********************************
  else 
  {
    fscanf(scn,"%d",&startyear);
    while((fscanf(scn,"%f",&t))!=EOF)
		{
      		if(t>999)
			endyear=(int)t;
    	}	
  }
  //**********************************
  fclose(scn);
  totalyears=endyear-startyear+1;
  
  if(e_flag) 
  {
    if(s_year==0) e_year=endyear;
    else if(e_year==0) e_year=s_year;
    s_year=startyear;
	    if(t_year>0 && (e_year-t_year+1)>startyear)
		{
	      s_year=e_year-t_year+1;
	    }
  }
  //**********************************
  else 
  {
    // Otherwise calculations are done from s_year.
    if(s_year==0) s_year=startyear;
    if(e_year==0) e_year=endyear;
	    if(t_year>0 && (s_year+t_year-1)<endyear)
		{
	      e_year=s_year+t_year-1;
	    }
  }
  //**********************************   
  if (setCalibrationStartYear == 1) 
  {
	 if (calibrationStartYear < startyear || calibrationStartYear > endyear) 
	 {
		calibrationStartYear = startyear;
	 }
  } 
  //**********************************
  else 
  {  /* no calibration start year set, so set it to startyear */
     calibrationStartYear = startyear;
  }
  currentCalibrationStartYear = calibrationStartYear;
  //**********************************
  if (setCalibrationEndYear == 1) 
  {
     if (calibrationEndYear < calibrationStartYear || calibrationEndYear > endyear) 
	 {
		calibrationEndYear = endyear;
	 }
  } 
  //**********************************
  else 
  {  /* no calibration end year set, so set it to startyear */
     calibrationEndYear = endyear;
  }
  
  currentCalibrationEndYear = calibrationEndYear;
  nStartYearsToSkip = currentCalibrationStartYear - startyear;
  nEndYearsToSkip = endyear - currentCalibrationEndYear;
  nCalibrationYears = currentCalibrationEndYear - currentCalibrationStartYear + 1;
  nStartPeriodsToSkip = nStartYearsToSkip * num_of_periods; 
  nEndPeriodsToSkip = nEndYearsToSkip * num_of_periods; 
  nCalibrationPeriods = nCalibrationYears * num_of_periods;
}// End set_flags
void pdsi::MonthlyPDSI() {
  int i;
  FILE *param;
  char filename[170];
  Monthly = true;
  SCMonthly = false;
  period_length = 1;
  num_of_periods = 12;
  currentCalibrationStartYear = startyear;
  currentCalibrationEndYear = endyear;
  nEndYearsToSkip = 0;
  nStartYearsToSkip = 0;
  nCalibrationYears = totalyears;
  nStartPeriodsToSkip = 0;
  nEndPeriodsToSkip = 0;
  nCalibrationPeriods = nCalibrationYears * num_of_periods;
  if(strlen(input_dir)>1)
    strcpy(filename,input_dir);
  else
    strcpy(filename,"./");
  strcat(filename,"parameter");
  
  if((param=fopen(filename,"r"))==NULL) 
  {
    exit(1);
  }
  GetParam(param);
  fclose(param);
  SumAll();
  // This outputs those sums to the screen
  for (i=0;i<num_of_periods;i++) 
  	{
    DEPSum[i] = ETSum[i] + RSum[i] - PESum[i] + ROSum[i];
    }
  
  DSSqr[i] = 0;
  CalcWBCoef();
  Calcd();
  CalcOrigK();
}
void pdsi::SCMonthlyPDSI() {
 int i;
  FILE *param;
  char filename[170];
  SCMonthly = true;
  Monthly = false; 
  if(initialize()<1){
    if(verbose > 0)
    if(verbose > 1){
    }
    return;
  }
  period_length = 1;
  num_of_periods = 12;
  currentCalibrationStartYear = calibrationStartYear;
  currentCalibrationEndYear = calibrationEndYear;
  nCalibrationYears = currentCalibrationEndYear - currentCalibrationStartYear + 1;
  nStartYearsToSkip = currentCalibrationStartYear - startyear;
  nEndYearsToSkip = endyear - currentCalibrationEndYear;
  nStartPeriodsToSkip = nStartYearsToSkip * num_of_periods; 
  nEndPeriodsToSkip = nEndYearsToSkip * num_of_periods; 
  nCalibrationPeriods = nCalibrationYears * num_of_periods;
  if(verbose > 1){   /* output calibration interval vars if max verbose is set */
  }
  if(strlen(input_dir)>1)
    strcpy(filename,input_dir);
  else
    strcpy(filename,"./");
  strcat(filename,"parameter");
  if((param=fopen(filename,"r"))==NULL) {
    exit(1);
  }
  GetParam(param);
  fclose(param);
  if(verbose>1)
  SumAll();
  // This outputs those sums to the screen
  if(verbose>1) {
  }
  for (i=0;i<num_of_periods;i++) {
    DEPSum[i] = ETSum[i] + RSum[i] - PESum[i] + ROSum[i];
    if(verbose>1) {
    }
    DSSqr[i] = 0;
  }
  CalcWBCoef();
  Calcd();
  CalcK();
  CalcZ();
  CalcDurFact(wetm, wetb, 1); 
  CalcDurFact(drym, dryb, -1); 
  if(verbose>1) { 
  }
  CalcX();
  Calibrate();
  if(verbose>1) {
	int i;
    for (i=0;i<num_of_periods;i++) {
    }
    for (i=0;i<num_of_periods;i++) {
      if (i==7) {
        number E, DE;
        E=SD/nCalibrationYears;
        DE=sqrt((SD2-E*SD)/(nCalibrationYears-1));
      }
    }
  }
}//end of SCMonthlyPDSI()
int pdsi::GetTemp(FILE *In, number *A, int max) {
  float t[52], t2[52], temp;
  int i, j, year, read, bad_weeks;
  char line[4096];
  char letter;
  for(i = 0; i < 52; i++)
    A[i] = 0;
  fgets(line,4096,In);
  read=sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",&year,&t[0],&t[1],&t[2],&t[3],&t[4],&t[5],&t[6],&t[7],&t[8],&t[9],&t[10],&t[11],&t[12],&t[13],&t[14],&t[15],&t[16],&t[17],&t[18],&t[19],&t[20],&t[21],&t[22],&t[23],&t[24],&t[25],&t[26],&t[27],&t[28],&t[29],&t[30],&t[31],&t[32],&t[33],&t[34],&t[35],&t[36],&t[37],&t[38],&t[39],&t[40],&t[41],&t[42],&t[43],&t[44],&t[45],&t[46],&t[47],&t[48],&t[49],&t[50],&t[51]);
  if(read == max+1)
  {
    //a full year's worth of data was read
    for(i = 0; i < max; i++)
	{
      t2[i] = t[i];
    }
  }
  //**********************************
  else
  {
    //check to see if it is the end of file
    if( (letter = fgetc(In)) != EOF ) 
	{
      for(i = 0 ; i < max - (read-1); i++)
		t2[i] = MISSING;
      for(i; i < max; i++)
		t2[i] = t[i - (max-read+1)];
      	ungetc(letter, In);
    }
    //**********************************
    else 
	{
      for(i = 0; i < read - 1; i++)
		t2[i] = t[i];
      for(i; i < max; i++)
		t2[i] = MISSING;
    }
  }
  //**********************************
  for(i = 0; i < num_of_periods; i++) 
  {
    bad_weeks = 0;
    temp = 0;
    for(j = 0; j < period_length; j++) 
	{
      if(t2[i*period_length + j] != MISSING)
		temp += t2[i*period_length + j];
      else
		bad_weeks++;
    }
    if(bad_weeks < period_length)
      A[i] = temp / (period_length - bad_weeks);
    else
      A[i] = MISSING;
  }
  //**********************************
  if(metric)
  {
	for(i = 0; i < num_of_periods; i++)
	{
	  if(A[i] != MISSING)
		A[i] = A[i] * (9.0/5.0) + 32;
	}
  }
  return year;
}
int pdsi::GetPrecip(FILE *In, number *A, int max) {
  float t[52], t2[52], temp;
  int i, j, year, read, bad_weeks;
  char line[4096];
  char letter;
  for(i = 0; i < 52; i++)
    A[i] = 0;
  fgets(line,4096,In);
  read=sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",&year,&t[0],&t[1],&t[2],&t[3],&t[4],&t[5],&t[6],&t[7],&t[8],&t[9],&t[10],&t[11],&t[12],&t[13],&t[14],&t[15],&t[16],&t[17],&t[18],&t[19],&t[20],&t[21],&t[22],&t[23],&t[24],&t[25],&t[26],&t[27],&t[28],&t[29],&t[30],&t[31],&t[32],&t[33],&t[34],&t[35],&t[36],&t[37],&t[38],&t[39],&t[40],&t[41],&t[42],&t[43],&t[44],&t[45],&t[46],&t[47],&t[48],&t[49],&t[50],&t[51]);
  if(read == max+1)
  {
    //a full year's worth of data was read
    for(i = 0; i < max; i++)
      t2[i] = t[i];
  }
  //**********************************
  else
  {
    //check to see if it is the end of file
    if( (letter = fgetc(In)) != EOF ) 
	{
      //it's not the end of the file
      //so place partial year's data at end of array
      for(i = 0 ; i < max - (read-1); i++)
		t2[i] = MISSING;
      for(i; i < max; i++)
		t2[i] = t[i - (max-read+1)];
      ungetc(letter, In);
    }
    else 
	{
      for(i = 0; i < read - 1; i++)
		t2[i] = t[i];
      for(i; i < max; i++)
		t2[i] = MISSING;
    }
  }
  //**********************************
  //now summaraize data in t2 into A
  for(i = 0; i < num_of_periods; i++) 
  {
    bad_weeks = 0;
    temp = 0;
    for(j = 0; j < period_length; j++) 
	{
      if(t2[i*period_length + j] != MISSING)
		temp += t2[i*period_length + j];
      else
		bad_weeks++;
    }
    //**********************************
    if(bad_weeks < period_length)
      A[i] = temp;
    else
      A[i] = MISSING;
  }
  //**********************************
  if(metric)
  {
	for(i = 0; i < num_of_periods; i++)
	{		
	  if(A[i] != MISSING)
		A[i] = A[i]/25.4;
	}
  }
  return year;
} 
void pdsi::GetParam(FILE * Param) {
  float scn1;
  fscanf(Param,"%f %f",&scn1);
  AWC=double(scn1);
  if(AWC>25.4){
  Ss = 25.4; 
  Su = AWC - Ss;
  }
  else 
  {	
    AWC =Ss;
    Su=0;
  }
}
void pdsi::CalcMonPE(int month, int year) {
  if (T[month] <= 92)
    PE=T[month];
}
void pdsi::Calcd() {
  FILE *fin;        // File pointer for the temp. input file potentials
  FILE *fout;       // File pointer for the temp. output file dvalue
  int per;           // The period in question
  int yr;           // The year in question
  number p;         // The precip for that period
  float scn1, scn2, scn3, scn4, scn5, scn6; // Temp. variables for fscanf
  char letter = ' ';
  int i = 0;
  number D_sum[52];
  number DSAct[52];
  number SPhat[52];
  for(i=0;i<52;i++)
  {
    D_sum[i] = 0.0;
	DSAct[i] = 0.0;
	SPhat[i] = 0.0;
  }
  // The potentials file is opened for reading in the previously stored values
  if((fin=fopen("potentials","r"))==NULL) { ; }
  if((fout=fopen("dvalue","w")) == NULL) { ; }
  //**********************************
  while(letter != '\n')
    letter = fgetc(fin);
  // This reads in the values previously stored in "potentials"
  //**********************************
  while(fscanf(fin,"%d %d %f %f %f %f %f %f", &yr, &per, &scn1, &scn2, &scn3, &scn4, &scn5, &scn6) != EOF) 
  {
    per = (per-1)/period_length;   //adjust the period # for use in arrays.
    p=scn1;
    PE=scn2;
    PR=scn3;
    PRO=scn4;
    PL=scn5;
	    if(p!=MISSING&&PE!=MISSING&&PR!=MISSING&&PRO!=MISSING&&PL!=MISSING)
		{
	      // Then the calculations for Phat and d are done
	      Phat=(Alpha[per]*PE)+(Beta[per]*PR)+(Gamma[per]*PRO)-(Delta[per]*PL);
	      d=p - Phat;
	      fprintf (fout, "%d %d %f\n", yr, (period_length*per)+1, d);
		     //**********************************
			 if (yr >= currentCalibrationStartYear && yr <= currentCalibrationEndYear) 
			   {
		        if(d < 0.0)
				D_sum[per] += -(d);
		        
				else
				D_sum[per] += d;
		        DSAct[per] += d;
		        DSSqr[per] += d*d;
		        SPhat[per] += Phat;
		       }
	    }
	    //**********************************
		else 
		{
	      d = MISSING;
	      fprintf (fout, "%d %d %f\n", yr, (period_length*per)+1, d);
	    }
  }
  //**********************************
  for(i = 0; i < num_of_periods; i++) 
  {
    D[i] = D_sum[i] / nCalibrationYears;
  }
  //**********************************
  // The files are closed 
  fclose(fin);
  fclose(fout);
}// end of Calcd()
void pdsi::CalcK() {
  number sums;        //used to calc k
  for(int per = 0; per < num_of_periods; per++)
  {
    if(PSum[per] + LSum[per] == 0)
      sums = 0;//prevent div by 0
    else
      sums = (PESum[per] + RSum[per] + ROSum[per]) / (PSum[per] + LSum[per]);
	
	//**********************************
    if(D[per] == 0)
      k[per] = 0.5;//prevent div by 0
    else
      k[per] = (1.5) * log10((sums + 2.8) / D[per]) + 0.5;
  } 
}//end of CalcK()
void pdsi::CalcOrigK() {
  int month, year;
  number sums;        //used to calc k
  float dtemp;
  DKSum = 0;
  FILE * inputd; // File pointer for input file dvalue
  // The dvalue file is open for reading. 
  if((inputd=fopen("dvalue","r"))==NULL) { ; } 
  //**********************************
  
  for(int per = 0; per < num_of_periods; per++)
  {
    if(PSum[per] + LSum[per] == 0)
      sums = 0;//prevent div by 0
    else
      sums = (PESum[per] + RSum[per] + ROSum[per]) / (PSum[per] + LSum[per]);
    if(D[per] == 0)
      k[per] = 0.5;//prevent div by 0
    else
      k[per] = (1.5) * log10((sums + 2.8) / D[per]) + 0.5;
	DKSum += D[per]*k[per];
  } 
  
  drym = .309;
  dryb = 2.691;
  wetm = drym;
  wetb = dryb;
  Prob = 0.0;
  X1 = 0.0;
  X2 = 0.0;
  X3 = 0.0;
  X = 0.0;
  V = 0.0;
  Q = 0.0;
  FILE* table;
  
  if((table=fopen("bigTablecalibrada.tbl","w"))==NULL){ ; }
		
  fprintf(table, "YEAR  MONTH     Z     %Prob     "); 
  fprintf(table, "X1       X2      X3\n"); 
  
  while((fscanf(inputd,"%d %d %f", &year, &month, &dtemp))!=EOF) 
  {
    PeriodList.insert(month);
    YearList.insert(year);
    d=dtemp;
    month--;
    K = (17.67/DKSum) * k[month];
    if(d != MISSING)
      Z = d*K;
    else
      Z = MISSING;
    ZIND.insert(Z);
    CalcOneX(table,month,year);
  }
  
  fclose(inputd);
  if(table)
  fclose(table);
  // Now that all calculations have been done they can be output to the screen
  if(verbose>1) 
  {
    int i;
    for (i=0;i<num_of_periods;i++) 
	{
      printf ("%4d %8.3f %8.2f %8.2f ", (period_length*i)+1, D[i], sqrt(DSSqr[i]/(totalyears-1)), DEPSum[i]/totalyears);
      if (i==7) 
	  {
        number E, DE;
        E=SD/totalyears;
        DE=sqrt((SD2-E*SD)/(totalyears-1));
      }
    }
  }
  
}//end of CalcOrigK()
void pdsi::CalcZ() {
  int year, per;
  float dtemp;
  llist tempZ,tempPer, tempyear;
  DKSum = 0.0; //sum of all D[i] and k[i]; used to calc K
  FILE * inputd; // File pointer for input file dvalue
  // The dvalue file is open for reading.
  if((inputd=fopen("dvalue","r"))==NULL) {
    if(verbose>0)
    exit(1);
  }
  for(per = 0; per < num_of_periods; per++)
	DKSum += D[per] * k[per];
  while((fscanf(inputd,"%d %d %f", &year, &per, &dtemp))!=EOF) { 
    PeriodList.insert(per);
    YearList.insert(year);
    per = (per-1)/period_length;//adjust for use in arrays
    d=dtemp;
    K = k[per];
    if(d != MISSING){
      Z = d*k[per];
    }
    else{
      Z = MISSING;
    } 
    ZIND.insert(Z);
  }
  fclose(inputd);
}//end of CalcZ()
void pdsi::Calibrate() {  
  llist tempZ, tempweek, tempyear;
  double cal_range;
  int size;
  llist tmpXlist; // temporary copy of Xlist, which is unscaled PDSI values
  //**********************************
  if ((setCalibrationStartYear == 1) || (setCalibrationEndYear == 1))
  {
     copy(tmpXlist,Xlist);  
	     for (int i=0; (i < nStartPeriodsToSkip) && (!tmpXlist.is_empty()) ; i++) 
		 {
	       tmpXlist.tail_remove(); /* remove periods before the start of the interval */
	     }
	     //**********************************
		 for (int i=0; (i < nEndPeriodsToSkip) && (!tmpXlist.is_empty()) ; i++) 
		 {
	       tmpXlist.head_remove(); /* remove periods before the end of the interval */
	     }
     size = tmpXlist.get_size();
	     if (nCalibrationPeriods != size) { ; }
     //calibrate using upper and lower 2% of values within the user-defined calibration interval
     cal_range = 4.0;
     dry_ratio = (-cal_range / tmpXlist.safe_percentile(0.02));
     wet_ratio = (cal_range / tmpXlist.safe_percentile(0.98));
/* SG 6/5/06: That ends the changes needed for implementing self-calibration intervals. */
  } 
  //**********************************
  else 
  {
     //calibrate using upper and lower 2% 
     cal_range = 4.0;
     dry_ratio = (-cal_range / Xlist.safe_percentile(0.02));
     wet_ratio = (cal_range / Xlist.safe_percentile(0.98));
  }
  //adjust the Z-index values
  //**********************************
  while(!ZIND.is_empty())
  {
    Z = ZIND.tail_remove();
    if(Z != MISSING)
	{
      if(Z >= 0)
		Z = Z * wet_ratio;
      else
		Z = Z * dry_ratio;
    }
    tempZ.insert(Z);
  }

  copy(ZIND,tempZ);
  CalcX();
}//end of calibrate()
void pdsi::CalcX() {
  llist tempZ,tempPer,tempYear;
  int year, per;
  FILE* table;
    table = fopen("bigTableautocalibrado.tbl","w");
    if(table == NULL){ ; }
    //**********************************
	else
	{ 
        fprintf(table, "YEAR  MONTH     Z     %Prob     ");  
        fprintf(table, "X1       X2      X3\n");  
    }  
  //**********************************
  while(!Xlist.is_empty())
    Xlist.head_remove();
  while(!XL1.is_empty())
    XL1.head_remove();
  while(!XL2.is_empty())
    XL2.head_remove();
  while(!XL3.is_empty())
    XL3.head_remove();
  while(!altX1.is_empty())
    altX1.head_remove();
  while(!altX2.is_empty())
    altX2.head_remove();
  while(!ProbL.is_empty())
    ProbL.head_remove();
  Prob = 0.0;
  X1 = 0.0;
  X2 = 0.0;
  X3 = 0.0;
  X = 0.0;
  V = 0.0;
  Q = 0.0;
  copy(tempZ, ZIND);
  copy(tempPer, PeriodList);
  copy(tempYear, YearList);
  
  while(!tempZ.is_empty())
  {	
    Z = tempZ.tail_remove();
    per = (int)tempPer.tail_remove();
    year = (int)tempYear.tail_remove();
    CalcOneX(table,per,year);
  }
  
  if(table)
  fclose(table);
}//end of CalcX()
void pdsi::CalcOneX(FILE* table, int period_number, int year) {
  number newV;    //These variables represent the values for 
  number newProb; //corresponding variables for the current period.
  number newX;    //They are kept seperate because many calculations
  number newX1;   //depend on last period's values.  
  number newX2;
  number newX3;
  number ZE;      //ZE is the Z value needed to end an established spell
  number m, b, c;
  flag wd;        //wd is a sign changing flag.  It allows for use of the same
  //**********************************
  if(X3>=0)
  {
    m = wetm;
    b = wetb;
  }
  //**********************************
  else
  {
    m = drym;
    b = dryb;
  }
  //**********************************
  c = 1 - (m / (m + b));
  if(Z != MISSING)
  {
  	//**********************************
  	//**********************************
    // This sets the wd flag by looking at X3
    if(X3>=0) wd=1;
    else wd=-1;
	    
		if(X3==0) 
		{
	      newX3=0;
	      newV=0;
	      newProb=0;
	      ChooseX(newX, newX1, newX2, newX3, bug);
	    }
	    //**********************************
	    else 
		{
		  newX3 = (c * X3 + Z/(m+b));
	      ZE = (m+b)*(wd*0.5 - c*X3);
	      Q=ZE+V;  
	      newV = Z - wd*(m*0.5) + wd*min(wd*V+tolerance,0);
	      //**********************************
		  if((wd*newV)>0) 
		  {
			newV=0;
			newProb=0;
			newX1=0;
			newX2=0;
			newX=newX3;
			//**********************************
				while(!altX1.is_empty())
				  altX1.head_remove();
				while(!altX2.is_empty())
				  altX2.head_remove();
	      }
	      else 
		  {
			newProb=(newV/Q)*100;
			if(newProb>=100-tolerance) 
			{
			  newX3=0;
			  newV=0;
			  newProb=100;
			}
			ChooseX(newX, newX1, newX2, newX3, bug);
	      }
	    }
	    
	//**********************************   
    if(table != NULL)
	{
      fprintf(table, "%5d %5d %7.2f %7.2f ",year,period_number,Z,newProb);
      fprintf(table, "%7.2f %7.2f %7.2f\n",newX1, newX2, newX3);
    }
    
	V = newV;
    Prob = newProb;
    X1 = newX1;
    X2 = newX2;
    X3 = newX3;
    Xlist.insert(newX);
    XL1.insert(X1);
    XL2.insert(X2);
    XL3.insert(X3);
    ProbL.insert(Prob);
  //**********************************
  //**********************************
  }
  
  else
  {
    if(table != NULL)
	{
      fprintf(table, "%5d %5d %7.2f %7.2f ",year,period_number,Z,MISSING);
      fprintf(table, "%7.2f %7.2f %7.2f\n",MISSING, MISSING, MISSING);
    }
    
	Xlist.insert(MISSING);
    XL1.insert(MISSING);
    XL2.insert(MISSING);
    XL3.insert(MISSING);
	ProbL.insert(MISSING);
  }
    
}//end of CalcOneX
void pdsi::SumAll() {
  FILE * fout;
  FILE * input_temp, *input_prec;
  char Temp[150], Precip[150];
  int actyear;
  number DEP=0;
  SD=0;
  SD2=0;
    /* SG 6/5/06: add variable to support a calibration interval */
  int nCalibrationPeriodsLeft = nCalibrationPeriods; /* init periods left */

  // Initializes the sums to 0;
  for(int i = 0; i < 52; i++) 
  {
    ETSum[i] = 0;
    RSum[i] = 0;
    LSum[i] = 0;
    ROSum[i] = 0;
    PSum[i] = 0;
    PESum[i] = 0;
    PRSum[i] = 0;
    PLSum[i] = 0;
    PROSum[i] = 0;
  }
  
  if((fout=fopen("potentials","w"))==NULL){ ; } 
  //**********************************
  else if(Monthly || SCMonthly)
  {

    if(strlen(input_dir)>1)
	{
      sprintf(Temp,"%s%s",input_dir,"monthly_T");
      sprintf(Precip,"%s%s",input_dir,"monthly_P");  
    }
    //**********************************
    else
	{
      strcpy(Temp,"monthly_T");
      strcpy(Precip,"monthly_P");
    }
    //**********************************
    if((input_temp=fopen(Temp,"r")) == NULL) { exit(1); }
    if((input_prec=fopen(Precip,"r")) == NULL) { exit(1); }
    //print column headers:
    fprintf(fout," Year  MONTH       P         PE         PR         PRO");
    fprintf(fout,"        PL        P-PE \n");
  }
  
  for(int year = 1; year <= totalyears; year++) 
  {

      actyear=GetTemp(input_temp, T, 12);
      GetPrecip(input_prec, P, 12);

    for(int per = 0; per < num_of_periods; per++) 
	{
      if(P[per] >= 0 && T[per] != MISSING)
	  {
		CalcMonPE(per,actyear);
		CalcPR();         // calculate Potential Recharge, Potential Runoff,
		CalcPRO();        // and Potential Loss
		CalcPL();
		CalcActual(per);  // Calculate Evapotranspiration, Recharge, Runoff,
		if(Weekly){
		  if (per > (17/period_length) && per < (35/period_length)) {
			DEP = DEP + P[per] + L - PE;
			if (per>(30/period_length) && per < (35/period_length)) {
			  SD=SD+DEP;
			  SD2=SD2+DEP*DEP;
			  DEP=0;
			}
		  }
		}
		else{
		  if (per > 4 && per < 8) {
			DEP = DEP + P[per] + L - PE;
			if (per == 7) {
			  SD=SD+DEP;
			  SD2=SD2+DEP*DEP;
			  DEP=0;
			}
		  }
		}
	if (year > nStartYearsToSkip) {  
	   if (nCalibrationPeriodsLeft > 0) { /* Within the calibration interval, so update Sums */
                // Reduce number of calibration years left by one for each year actually summed 
	        nCalibrationPeriodsLeft--; 
		ETSum[per] += ET;
		RSum[per] += R;
		ROSum[per] += RO;
		LSum[per] += L;
		PSum[per] += P[per];
		PESum[per] += PE;
		PRSum[per] += PR;
		PROSum[per] += PRO;
		PLSum[per] += PL;
	   }
	}
		fprintf(fout,"%5d %5d %10.6f ",actyear,(period_length*per)+1,P[per]);
		fprintf(fout,"%10.6f %10.6f %10.6f %10.6f ",PE,PR,PRO,PL);
		fprintf(fout,"%10.6f\n",P[per]-PE);
      }//matches if(P[per]>= 0 && T[per] != MISSING)
      else {
		fprintf(fout,"%5d %5d %f ",actyear, (period_length*per)+1,MISSING);
		fprintf(fout,"%10.6f %10.6f %10.6f ", MISSING, MISSING, MISSING);
		fprintf(fout,"%10.6f %10.6f\n",MISSING, MISSING);
      }
    }//end of period loop
  }//end of year loop
  fclose(fout);
  fclose(input_temp);
  fclose(input_prec);
}
void pdsi::CalcPR() {
  PR = AWC - (Su + Ss);
}
void pdsi::CalcPRO() {
  PRO = Ss + Su;
}
void pdsi::CalcPL() {
  if(Ss >= PE)
    PL = PE;
  else {
    PL = ((PE - Ss) * Su) / (AWC) + Ss;
    if(PL > PRO)  // If PL>PRO then PL>water in the soil.  This isn't
      PL = PRO;   // possible so PL is set to the water in the soil
  }
}
void pdsi::CalcActual(int per) {
  number R_surface = 0.0;   // recharge of the surface layer
  number R_under = 0.0;    // recharge of the underlying layer
  number surface_L = 0.0;   // loss from surface layer
  number under_L = 0.0;    // loss from underlying layer
  number new_Su, new_Ss;    // new soil moisture values
  if(P[per] >= PE) {
    ET = PE;   // Enough moisture for all potential evapotranspiration to occur
    L = 0.0;   // with no actual loss of soil moisture
    if((P[per] - PE) > (1.0 - Ss)) {
      R_surface = 1.0 - Ss;
      new_Ss = 1.0;
      if((P[per] - PE - R_surface) < ((AWC - 1.0) - Su)) {
		R_under = (P[per] - PE - R_surface);
		RO = 0.0;
      }
      else {
		R_under = (AWC - 1.0) - Su;
		RO = P[per] - PE - (R_surface + R_under);
      }
      new_Su = Su + R_under;
      R = R_surface + R_under;//total recharge
    }
    else {
      // There is only enough moisture to recharge some of the top layer.
      R = P[per] - PE;
      new_Ss = Ss + R;
      new_Su = Su;
      RO = 0.0;
    }
  }// End of if(P[per] >= PE)
  else {
    // The evapotranspiration is greater than the precipitation received.  This
    // means some moisture loss will occur from the soil.
    if(Ss > (PE - P[per])) {
      surface_L = PE - P[per];
      under_L = 0.0;
      new_Ss = Ss - surface_L;
      new_Su = Su;
    }
    else {
      // The top layer is drained, so the underlying layer loses moisture also.
      surface_L = Ss;
      under_L = (PE - P[per] - surface_L) * Su / AWC;
      if(Su < under_L)
		under_L = Su;
      new_Ss = 0.0;
      new_Su = Su - under_L;
    }
    R = 0;// No recharge occurs
    L = under_L + surface_L;// Total loss
    RO = 0.0;// No extra moisture so no runoff
    ET = P[per] + L;// Total evapotranspiration
  }
  Ss = new_Ss;//update soil moisture values
  Su = new_Su;
  
}//end of CalcActual(int per)
void pdsi::CalcWBCoef() {
  FILE *wb;
  // The coefficients are calculated by per
  for (int per=0; per < num_of_periods; per++) {  
    //calculate alpha:
    if(PESum[per] != 0.0)
      Alpha[per] = ETSum[per] / PESum[per];
    else if(ETSum[per] == 0.0)
      Alpha[per] = 1.0;
    else
      Alpha[per] = 0.0;
    //calculate beta:
    if(PRSum[per] != 0.0)
      Beta[per] = RSum[per] / PRSum[per];
    else if(RSum[per] == 0.0)
      Beta[per] = 1.0;
    else
      Beta[per] = 0.0;
    //calculate gamma:
    if(PROSum[per] != 0.0)
      Gamma[per] = ROSum[per] / PROSum[per];
    else if(ROSum[per] == 0.0)
      Gamma[per] = 1.0;
    else
      Gamma[per] = 0.0;
    //calculate delta:
    if(PLSum[per] != 0.0)
      Delta[per] = LSum[per] / PLSum[per];
    else 
      Delta[per] = 0.0;
  }
 //output water balance coefficients
    if((wb=fopen("WB.tbl","w"))==NULL){
    fprintf(wb, "ETPRIOD   ALPHA     BETA    GAMMA    DELTA\n");			
										}
	if(verbose>1)
      printf ("\nETPRIOD   ALPHA     BETA    GAMMA    DELTA\n");
    for(int i = 0; i < num_of_periods; i++){
      fprintf(wb, "%3d %10.4f %8.4f %8.4f %8.4f \n", (period_length*i)+1, Alpha[i], Beta[i], Gamma[i], Delta[i]);
      if(verbose>1)
		printf ("%3d %10.4f %8.4f %8.4f %8.4f \n", (period_length*i)+1, Alpha[i], Beta[i], Gamma[i], Delta[i]);
    }
    fclose(wb);
}//end CalcWBCoef()
void pdsi::Write() {
  char full_path[128];
  if(Weekly){
    sprintf(full_path, "weekly/%d/",period_length);
    Write(full_path);
  }
  else if(Monthly)
    Write("monthly/original");
  else if(SCMonthly)
    Write("monthly/self_cal");
  
}
void pdsi::Write(char *directory) {
  unsigned int i = 0;
  int e = 0;                            //error flag
  char base[128];                       //string for base dir path
  char my_dir[128];                     //string for directory
  char full_path[256];                  //string for full path to dir
  strcpy(my_dir,directory);
  //make sure last character of directory is '/'
  if(my_dir[strlen(my_dir)-1] != '/')
    strcat(my_dir,"/");
  if(strlen(output_dir)>1)
    strcpy(base,output_dir);
  else
    strcpy(base,"./");
  sprintf(full_path,"%s%s",base,my_dir);
  switch(create_dir(full_path)){
  case -1:
    if(verbose > 0)
      printf("Error creating directory: %s \n",full_path);
    for(i = 0; i < strlen(my_dir); i++){
      if(my_dir[i] == '/')
		my_dir[i] = '_';
    }
    sprintf(full_path,"%s%s",base,my_dir);
    if(verbose > 0)
      printf("Output files will now have the prefix: %s \n\n",my_dir);
    break;
  case 0:
    //the directory was successfully created
    if(verbose > 1)
      printf("Outputting to the new directory: %s\n\n",full_path);
    break;
  case 1:
    //the directory already exists
    if(verbose > 1)
      printf("Outputting to existing directory: %s\n\n",full_path);
    break;
  }

    FinalWrite(full_path);
    MoveFiles(full_path);
  
}
void pdsi::FinalWrite(char* dir) {
  int cyear=startyear;
  int prev, cur, saved_per, change_per;
  number x, x1, x2, x3, p, wp, ph, z;
  llist tempX, tempX1, tempX2, tempX3, tempZ, tempP,tempweek;
  FILE *pdsi0, *pdsi1, *phdi0, *phdi1, *zind0, *zind1, *wplm0, *wplm1;
  char filename[80];
  if(Xlist.is_empty()){
    if(verbose > 1)
      printf("No PDSI values have been calculated.\n");
    return;
  }
  //remove any existing output files:
  sprintf(filename,"%s%s",dir,"PDSI.tbl");
  remove(filename);
  sprintf(filename,"%s%s",dir,"PDSI.clm");
  remove(filename);
  sprintf(filename,"%s%s",dir,"PHDI.tbl");
  remove(filename);
  sprintf(filename,"%s%s",dir,"PHDI.clm");
  remove(filename);
  sprintf(filename,"%s%s",dir,"WPLM.tbl");
  remove(filename);
  sprintf(filename,"%s%s",dir,"WPLM.clm");
  remove(filename);
  sprintf(filename,"%s%s",dir,"ZIND.tbl");
  remove(filename);
  sprintf(filename,"%s%s",dir,"ZIND.clm");
  remove(filename);
  if(output_mode==0 || output_mode==2) {
    //open files
    sprintf(filename,"%s%s",dir,"PDSI.tbl");
    pdsi0=fopen(filename,"w");
    sprintf(filename,"%s%s",dir,"PHDI.tbl");
    phdi0=fopen(filename,"w");
    sprintf(filename,"%s%s",dir,"WPLM.tbl");
    wplm0=fopen(filename,"w"); 
    sprintf(filename,"%s%s",dir,"ZIND.tbl");
    zind0=fopen(filename,"w");
  }
  // If the output_mode is 1 or 2 the column files are opened for writing
  if(output_mode==1 || output_mode==2) {
    //open files
    sprintf(filename,"%s%s",dir,"PDSI.clm"); 
    pdsi1=fopen(filename,"w"); 
    sprintf(filename,"%s%s",dir,"PHDI.clm"); 
    phdi1=fopen(filename,"w");
    sprintf(filename,"%s%s",dir,"WPLM.clm");
    wplm1=fopen(filename,"w");
    sprintf(filename,"%s%s",dir,"ZIND.clm");
    zind1=fopen(filename,"w");
  }
  copy(tempX,Xlist);
  copy(tempX1,XL1);
  copy(tempX2,XL2);
  copy(tempX3,XL3);
  copy(tempZ,ZIND);
  copy(tempP,ProbL);
  copy(tempweek, PeriodList);
  change_per = 1;
  prev = 0;
  // The temp lists are then emptied and output
  while(!tempX.is_empty() && cyear<=e_year) {
    if(change_per)
      cur = (int)tempweek.tail_remove();
    else
      cur = saved_per;
    if((prev == 0 && cur == 1) || cur == prev + period_length || cur == 1 && ( (Weekly && prev == (52-period_length+1)) || ((Monthly||SCMonthly) && prev == (12-period_length+1)) ) ) {
      change_per = 1;
      x=tempX.tail_remove();
      x1=tempX1.tail_remove();
      x2=tempX2.tail_remove();
      x3=tempX3.tail_remove();
      p=tempP.tail_remove();
      z=tempZ.tail_remove();
      p=p/100;
      ph=x3;
      if (x3==0) {
		ph=x;
		wp=x1;
		if (-x2>(x1+tolerance))
		  wp=x2;
      }
      else if (p>(0+tolerance/100) && p<(1-tolerance/100)) {
		if (x3 < 0) 
		  // X3 is negative so WPLM is weighted average of X3 and X1
		  wp=(1-p)*x3 + p*x1;
		else
		  // X3 is positive so WPLM is weighted average of X3 and X2
		  wp=(1-p)*x3 + p*x2;
      }
      else
		wp=x3;
    }
    else{
      if(verbose>1)
		printf("missing data before period %d in %d\n",cur,cyear);
      saved_per = cur;
      cur = prev + period_length; 
      change_per = 0;
      x = MISSING;
      ph = MISSING;
      z = MISSING;
      wp = MISSING;
    }
    if(cur==1){
      if((output_mode==0 || output_mode==2) && cyear>=s_year){
		fprintf (pdsi0,"%5d",cyear);
		fprintf (phdi0,"%5d",cyear);
		fprintf (zind0,"%5d",cyear);
		fprintf (wplm0,"%5d",cyear);
      }
    }
    // Outputs the current values in a column format: Year Week Value
    if((output_mode==1 || output_mode==2) && cyear>=s_year) {
      fprintf (pdsi1,"%5d%5d%7.2f\n", cyear, cur, x);
      fprintf (phdi1,"%5d%5d%7.2f\n", cyear, cur, ph);
      fprintf (zind1,"%5d%5d%7.2f\n", cyear, cur, z);
      fprintf (wplm1,"%5d%5d%7.2f\n", cyear, cur, wp);
    }
    if((output_mode==0 || output_mode==2) && cyear>=s_year) {
      fprintf (pdsi0,"%7.2f", x);
      fprintf (phdi0,"%7.2f", ph);
      fprintf (zind0,"%7.2f", z);
      fprintf (wplm0,"%7.2f", wp);
    }
    if(Weekly && cur>=52-period_length+1 || (Monthly||SCMonthly) && cur>=12-period_length+1){
      // Outputs a newline after a years worth of data in the table files
      if((output_mode==0 || output_mode==2) && cyear>=s_year) {
		fprintf (pdsi0,"\n");
		fprintf (phdi0,"\n");
		fprintf (zind0,"\n");
		fprintf (wplm0,"\n");
      }
      cyear++;
      cur = 1 - period_length;
    }
    prev = cur;
  }
  if(verbose > 0){
    if(Weekly){
      if(period_length == 1)
		printf("%4s Self-Calibrated Weekly PDSI written to %s\n","*",dir);
      else if(period_length == 2)
		printf("%4s Self-Calibrated 2-Week PDSI written to %s\n","*",dir);
      else if(period_length == 4)
		printf("%4s Self-Calibrated 4-Week PDSI written to %s\n","*",dir);
      else if(period_length == 13)
		printf("%4s Self-Calibrated 13-Week PDSI written to %s\n","*",dir);
    }
    if(Monthly)
      printf("%4s Monthly PDSI written to %s\n","*",dir);
    if(SCMonthly)
      printf("%4s Self-Calibrated Monthly PDSI written to %s\n","*",dir);
  }
}// void pdsi::FinalWrite()
void pdsi::MoveFiles(char* dir) {
  char filename[170];
  //remove any files that might already be in the director
  sprintf(filename,"%s%s",dir,"potentials");
  remove(filename);
  sprintf(filename,"%s%s",dir,"dvalue");
  remove(filename);
  sprintf(filename,"%s%s",dir,"bigTable.tbl");
  remove(filename);
  sprintf(filename,"%s%s",dir,"WB.tbl");
  remove(filename);
  //move files
  if(extra == 3 || extra == 9){
    //move potentials and dvalue
    sprintf(filename,"%s%s",dir,"potentials");
    if(rename("potentials",filename) == -1){
      if(ReprintFile("potentials",filename) == -1 && verbose > 0)
		printf("Unable to rename potentials file as \"%s\".\n",filename);
    }
    sprintf(filename,"%s%s",dir,"dvalue"); 
    if(rename("dvalue",filename) == -1) { 
      if(ReprintFile("dvalue",filename) == -1 && verbose > 0) 
		printf("Unable to rename dvalue file as \"%s\".\n",filename); 
    } 
  }
  if(extra == 2 || extra == 9){
    //move bigTable.tbl
    sprintf(filename,"%s%s",dir,"bigTable.tbl"); 
    if(rename("bigTable.tbl",filename) == -1){ 
      if(ReprintFile("bigTable.tbl",filename) == -1 && verbose > 0) 
        printf("Unable to rename bigTable.tbl file as \"%s\".\n",filename); 
    } 
  }  
  if(extra == 1 || extra == 9){
    //move WB.tbl
    sprintf(filename,"%s%s",dir,"WB.tbl"); 
    if(rename("WB.tbl",filename) == -1 ) { 
      if(ReprintFile("WB.tbl",filename) == -1 && verbose > 0) 
        printf("Unable to rename WB.tbl file as \"%s\".\n",filename); 
    } 
  }
}
int pdsi::ReprintFile(char* src, char *des){
  FILE *in, *out;
  char line[4096];
  in = fopen(src,"r");
  if(in == NULL)
    return -1;
  out = fopen (des,"w");
  if(out == NULL)
    return 1;
  while(fgets(line, 4096, in))
    fprintf(out, "%s",line);
  if(feof(in)){
    fclose(out);
    fclose(in);
    remove(src);
    return 1;
  }
  else {
    fclose(in);
    fclose(out);
    return -1;
  }
}
void pdsi::Backtrack(number X1, number X2) {
  number num1,num2;
  node * ptr=NULL;
  num1=X1;
  while (!altX1.is_empty() && !altX2.is_empty()) {
    if (num1>0) {
      num1=altX1.head_remove();
      num2=altX2.head_remove();
    }
    else {
      num1=altX2.head_remove();
      num2=altX1.head_remove();
    }
    if (-tolerance<=num1 && num1<=tolerance) num1=num2;
    ptr=Xlist.set_node(ptr,num1);
  }
}//end of backtrack()
void pdsi::ChooseX(number& newX, number& newX1, number& newX2, number& newX3, int bug){
  number m, b;
  number wetc, dryc;

  if(X3>=0){
    m = wetm;
    b = wetb;
  }
  else{
    m = drym;
    b = dryb;
  }
  wetc = 1 - (wetm / (wetm+wetb));
  dryc = 1 - (drym / (drym+wetb));
  newX1 = (wetc*X1 + Z/(wetm+wetb));
  if(newX1 < 0)
    newX1 = 0;
  newX2 = X2;
  if(bug==0){
    newX2 = (dryc*X2 + Z/(drym+dryb));
    if(newX2 > 0)
      newX2 = 0;
  }
  if((newX1 >= 0.5)&&(newX3 == 0)){
    Backtrack(newX1, newX2);
    newX = newX1;
    newX3 = newX1;
    newX1 = 0;
  }
  else{
    newX2 = (dryc*X2 + Z/(drym+dryb));
    if(newX2 > 0)
      newX2 = 0;
    
    if((newX2 <= -0.5)&&(newX3 == 0)){
      Backtrack(newX2, newX1);
      newX = newX2;
      newX3 = newX2;
      newX2 = 0;
    }
    else if(newX3 == 0) {
      if(newX1 == 0){
		Backtrack(newX2, newX1);
		newX = newX2;
      }
      else if(newX2 == 0){
		Backtrack(newX1, newX2);
		newX = newX1;
      }
      else{
		altX1.insert(newX1);
		altX2.insert(newX2);
		newX = newX3;
      }
    }

    else{
      //store X1 and X2 in their linked lists for possible use later
      altX1.insert(newX1);
      altX2.insert(newX2);
      newX = newX3;
    }
  }
}//end of chooseX
void pdsi::CalcDurFact(number &slope, number &intercept, int sign){
  int num_list = 10;
  number sum[10];
  int length[10];
  int i;
  if(Weekly){
    if(period_length==1){
      length[0]=13;
      length[1]=26;
      length[2]=39;
      length[3]=52;
      length[4]=78;
      length[5]=104;
      length[6]=130;
      length[7]=156;
      length[8]=182;
      length[9]=208;
    }
    else if(period_length==2){
      length[0]=6;
      length[1]=13;
      length[2]=19;
      length[3]=26;
      length[4]=39;
      length[5]=52;
      length[6]=65;
      length[7]=78;
      length[8]=91;
      length[9]=104;
    }
    else if(period_length==4){
      length[0]=3;
      length[1]=6;
      length[2]=10;
      length[3]=13;
      length[4]=20;
      length[5]=26;
      length[6]=33;
      length[7]=39;
      length[8]=46;
      length[9]=52;
    }
    else if(period_length==13){
      length[0]=2;
      length[1]=3;
      length[2]=4;
      length[3]=5;
      length[4]=6;
      length[5]=8;
      length[6]=10;
      length[7]=12;
      length[8]=14;
      length[9]=16;
    }
  }
  else{
    length[0]=3;
    length[1]=6;
    length[2]=9;
    length[3]=12;
    length[4]=18;
    length[5]=24;
    length[6]=30;
    length[7]=36;
    length[8]=42;
    length[9]=48;
  }
  for(i = 0; i < num_list; i++){
    sum[i] = get_Z_sum(length[i],sign);
    //printf("sum[%d] = %f\n",i,sum[i]);
  }
  if(verbose > 1){
    printf("Points used in linear regression for Duration Factors:\n");
    for(i=0;i<num_list;i++)
      printf("%7d  ",length[i]);
    printf("\n");
    for(i=0;i<num_list;i++)
      printf("%7.2f  ",sum[i]);
    printf("\n");
  }
  LeastSquares(length, sum, num_list, sign, slope, intercept);
  slope = slope / (sign*4);
  intercept = intercept / (sign*4);
}//end of CalcDurFact()
number pdsi::get_Z_sum(int length, int sign) {
  number sum, max_sum, z;
  llist tempZ, list_to_sum, list_of_sums;

  number highest_reasonable;
  number percentile;
  number reasonable_tol = 1.25;
/* SG 6/5/06: Add variable to implement user-defined calibration interval */
  int nCalibrationPeriodsLeft;
  copy(tempZ,ZIND); 
  sum = 0;
  for (int i=0; (i < nStartPeriodsToSkip) && (!tempZ.is_empty()) ; i++) {
    tempZ.tail_remove(); /* remove periods before the start of the interval */
  }
   nCalibrationPeriodsLeft = nCalibrationPeriods; /* init periods left */
  //first fill the list to be summed  
  for(int i = 0; i < length; i++){
    if(tempZ.is_empty()){
      printf("Error: tempZ is empty.\n");
      i = length;
    }
    else {
      z = tempZ.tail_remove();
      nCalibrationPeriodsLeft--; /* reduce by one period for each remove */
      if(z != MISSING){
		sum += z;
		list_to_sum.insert(z);
      }
      else{
		i--;
      }
    }
    //printf("i = %d  z= %f  sum = %f\n",i,z,sum);
  }
  max_sum = sum;
  list_of_sums.insert(sum);
  while(!tempZ.is_empty() && nCalibrationPeriodsLeft > 0){
    z = tempZ.tail_remove();
    nCalibrationPeriodsLeft--; /* reduce by one period for each remove */
    if(z != MISSING){
      sum -= list_to_sum.tail_remove();
      sum += z;
      list_to_sum.insert(z);
      list_of_sums.insert(sum);
    }
    if(sign * sum > sign * max_sum)
      max_sum = sum;
  }
  highest_reasonable = 0; 
  if(sign == 1)
    percentile = list_of_sums.safe_percentile(.98);
  if(sign == -1)
    percentile = list_of_sums.safe_percentile(.02);
  while(!list_of_sums.is_empty()){
    sum = list_of_sums.tail_remove();
    if(sign * sum > 0 ){
      if( (sum / percentile) < reasonable_tol ) {
		if(sign * sum > sign * highest_reasonable )
		  highest_reasonable = sum;
      }
    }
  }
  if(sign == -1)
    return max_sum;
  else if(sign == 1)
    //return max_sum;
    return highest_reasonable;
  else
    return MISSING;
}//end of get_Z_sum()
void pdsi::LeastSquares(int *x, number *y, int n, int sign, number &slope, number &intercept) {
  number sumX, sumX2, sumY, sumY2, sumXY;
  number SSX, SSY, SSXY;
  number xbar, ybar;
  number correlation = 0;
  number c_tol = 0.85;
  number max = 0;
  number max_diff = 0;
  int max_i = 0;
  number this_x, this_y;
  int i;
  sumX = 0; sumY = 0; sumX2 = 0; sumY2 = 0; sumXY = 0;
  for(i = 0; i < n; i++){
    this_x = x[i];
    this_y = y[i];
    sumX += this_x;
    sumY += this_y;
    sumX2 += this_x * this_x;
    sumY2 += this_y * this_y;
    sumXY += this_x * this_y;
  }
  xbar = sumX/n;
  ybar = sumY/n;
  SSX = sumX2 - (sumX * sumX)/n;
  SSY = sumY2 - (sumY * sumY)/n;
  SSXY = sumXY - (sumX * sumY)/n;
  correlation = SSXY / (sqrt(SSX) * sqrt(SSY));
  if(verbose > 1 && (sign*correlation) < c_tol ){
    printf("original correlation = %.4f \n",correlation);
  }
  i = n - 1;
  while((sign*correlation) < c_tol && i > 3){
    this_x = x[i];
    this_y = y[i];

    sumX -= this_x;
    sumY -= this_y;
    sumX2 -= this_x * this_x;
    sumY2 -= this_y * this_y;
    sumXY -= this_x * this_y;

    SSX = sumX2 - (sumX * sumX)/i;
    SSY = sumY2 - (sumY * sumY)/i;
    SSXY = sumXY - (sumX * sumY)/i;

    xbar = sumX/i;
    ybar = sumY/i;

    correlation = SSXY / (sqrt(SSX) * sqrt(SSY));
    i--;
  }
  if(verbose > 1){
    printf("final correlation =  %.4f\n\n",correlation);
  }
  slope = SSXY / SSX;
  n = i+1;
  i = 0;
  for(i; i < n; i++){
    if(sign*(y[i] - slope * x[i]) > sign*max_diff){
      max_diff = y[i] - slope * x[i];
      max_i = i;
      max = y[i];
    }
  }
  intercept = max - slope*x[max_i];
}//end of LeastSquares()
number pdsi::getPDSI(int period, int year) {
  return getValue(Xlist, period, year);
}
number pdsi::getZIND(int period, int year) {
  return getValue(ZIND, period, year);
}
number pdsi::getWPLM(int period, int year) {
  number x1,x2,x3,p,wp;
  x1 = getValue(XL1, period, year);
  x2 = getValue(XL2, period, year);
  x3 = getValue(XL3, period, year);
  p = getValue(ProbL, period, year);
  if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING)  
	wp = MISSING;  
  else{  
	p = p / 100;
	if (x3==0) {
	  wp=x1;
	  if (-x2>(x1+tolerance))
		wp=x2;
	}
	else if (p>(0+tolerance/100) && p<(1-tolerance/100)) {
	  if (x3 < 0) 
		// X3 is negative so WPLM is weighted average of X3 and X1
		wp=(1-p)*x3 + p*x1;
	  else
		// X3 is positive so WPLM is weighted average of X3 and X2
		wp=(1-p)*x3 + p*x2;
	}
	else
	  wp=x3;
  }
  return wp;
}
number pdsi::getPHDI(int period, int year) {
  number x, x3;
  x = getValue(Xlist, period, year);
  x3 = getValue(XL3, period, year);
  if(x == MISSING || x3 == MISSING) 
	return MISSING;
  if(x3==0)
	return x;
  else
	return x3;
}
number pdsi::getValue(llist &List, int period, int year) {
  llist tempPer, tempYear, tempList;
  number per, yr, val;
  bool loop_exit = false;
  copy(tempList, List);
  copy(tempPer, PeriodList);
  copy(tempYear, YearList);
  while(! loop_exit) {
	if(tempList.is_empty())
	  loop_exit = true;
	if(tempPer.is_empty())
	  loop_exit = true;
	if(YearList.is_empty())
	  loop_exit = true;
	val = tempList.head_remove();
	per = tempPer.head_remove();
	yr = tempYear.head_remove();

	if(yr == year && per == period)
	  return val;
  }
  return MISSING;
}
number* pdsi::getYearArray(int &size) {
  size = YearList.get_size();
  return YearList.returnArray();
}
number* pdsi::getPerArray(int &size) {
  size = PeriodList.get_size();
  return PeriodList.returnArray();
}
number* pdsi::getPDSIArray(int &size) {
  size = Xlist.get_size();
  return Xlist.returnArray();
}
number* pdsi::getPDSIArray(int start_per, int start_yr, int end_per, int end_yr, int &size) {
  return getSubArray(Xlist, start_per, start_yr, end_per, end_yr, size);
}
number* pdsi::getZINDArray(int &size){
  size = ZIND.get_size();
  return ZIND.returnArray();
}
number* pdsi::getZINDArray(int start_per, int start_yr, int end_per, int end_yr, int &size) {
  return getSubArray(ZIND, start_per, start_yr, end_per, end_yr, size);
}
number* pdsi::getPHDIArray(int &size) {
  number *x, *x3;
  x = Xlist.returnArray();
  x3 = XL3.returnArray();
  size = Xlist.get_size();
  number *A = new number[size];
  if(A == NULL){
	size = 0; 
	return A;
  }
  for(int i = 0; i < size; i++){
	if(x[i] != MISSING){
	  if(x3[i]==0)
		A[i] = x[i];
	  else
		A[i] = x3[i];
	}
	else
	  A[i] = MISSING;
  }
  delete [] x;
  delete [] x3;
  return A;
}
number* pdsi::getPHDIArray(int start_per, int start_yr, int end_per, int end_yr, int &size) {
  number *x, *x3; 
  int tempsize;
  x = getSubArray(Xlist, start_per, start_yr, end_per, end_yr, size);
  x3 = getSubArray(XL3, start_per, start_yr, end_per, end_yr, tempsize); 
  number *A = new number[size];
  if(A == NULL){
	size = 0;
	return A;
  }
  for(int i = 0; i < size; i++){ 
    if(x[i] != MISSING){ 
      if(x3[i]==0) 
        A[i] = x[i]; 
      else 
        A[i] = x3[i]; 
    } 
    else 
      A[i] = MISSING; 
  }
  delete [] x;
  delete [] x3;
  return A;
}
number* pdsi::getWPLMArray(int &size) {
  number *A;
  number *x1Array,*x2Array,*x3Array,*pArray; 
  number x1, x2, x3, p, wp;
  x1Array = XL1.returnArray();
  x2Array = XL2.returnArray();
  x3Array = XL3.returnArray();
  pArray = ProbL.returnArray();
  size = XL1.get_size();
  
  A = new number[size];
  if(A == NULL){
	size = 0;
	return A;
  }

  for(int i = 0; i < size; i++){
	x1 = x1Array[i];
	x2 = x2Array[i];
	x3 = x3Array[i];
	p = pArray[i];

	if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING)
	  wp = MISSING;
	else{
	  p = p / 100; 
	  if (x3==0) { 
		wp=x1; 
		if (-x2>(x1+tolerance)) 
		  wp=x2; 
	  } 
	  else if (p>(0+tolerance/100) && p<(1-tolerance/100)) { 
		if (x3 < 0)  
		  // X3 is negative so WPLM is weighted average of X3 and X1 
		  wp=(1-p)*x3 + p*x1; 
		else 
		  // X3 is positive so WPLM is weighted average of X3 and X2 
		  wp=(1-p)*x3 + p*x2; 
	  } 
	  else 
		wp=x3; 
	}
	A[i] = wp;
  }
  delete [] x1Array;
  delete [] x2Array;
  delete [] x3Array;
  delete [] pArray;
  return A;
}
number* pdsi::getWPLMArray(int start_per, int start_yr, int end_per, int end_yr, int &size) {
  number *A; 
  number *x1Array,*x2Array,*x3Array,*pArray;  
  number x1, x2, x3, p, wp; 
  int tempsize;
  x1Array = getSubArray(XL1, start_per, start_yr, end_per, end_yr, size);
  x2Array = getSubArray(XL2, start_per, start_yr, end_per, end_yr, tempsize);
  x3Array = getSubArray(XL3, start_per, start_yr, end_per, end_yr, tempsize);
  pArray = getSubArray(ProbL, start_per, start_yr, end_per, end_yr, tempsize);
  A = new number[size]; 
  if(A == NULL){ 
    size = 0; 
    return A; 
  } 
 
  for(int i = 0; i < size; i++){ 
    x1 = x1Array[i]; 
    x2 = x2Array[i]; 
    x3 = x3Array[i]; 
    p = pArray[i];
 
	if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING) 
      wp = MISSING; 
	else{ 
	  p = p / 100; 
	  if (x3==0) {  
		wp=x1;  
		if (-x2>(x1+tolerance))  
		  wp=x2;  
	  }  
	  else if (p>(0+tolerance/100) && p<(1-tolerance/100)) {  
		if (x3 < 0)   
		  wp=(1-p)*x3 + p*x1;  
		else  
		  wp=(1-p)*x3 + p*x2;  
	  }  
	  else  
		// There is an established spell without possibility of end meaning the  
		// WPLM is simply X3  
		wp=x3;  
	}
    A[i] = wp; 
  } 
  delete [] x1Array;
  delete [] x2Array;
  delete [] x3Array;
  delete [] pArray;
  return A; 
}
number* pdsi::getSubArray(llist &List, int start_per, int start_yr, int end_per, int end_yr, int &size) {
  llist temp;
  number *Array, *year, *period;
  int i,j;
  int cur_per, cur_yr;
  int per_len=0;
  int num_missing;
  Array = List.returnArray();
  year = YearList.returnArray();
  period = PeriodList.returnArray();
  for(j = 0; j < PeriodList.get_size(); j++){
    if(period[j] > per_len)
      per_len = (int)period[j];
  }   
  if( (start_yr > year[0]) ||
      ( (start_yr == year[0]) && (start_per > period[0]) ) ) { 
    i = 0;
    while( ( ( year[i] < start_yr ) || 
		   ( year[i] == start_yr && period[i] < start_per ) ) &&
		 ( i < List.get_size() ) ) {
	i++;
    }
    while( ( ( year[i] < end_yr ) ||
		   ( year[i] == end_yr && period[i] <= end_per ) ) &&
		 ( i < List.get_size() ) ) {
	temp.insert(Array[i]);
	i++;
    }
    if(i == List.get_size()){
      cur_yr = (int)year[i-1];
      cur_per = (int)period[i-1];
      if((cur_per%per_len) == 0){
        cur_per = 1;
        cur_yr++;
      }
      else
        cur_per++;   
      while( (cur_yr < end_yr) ||
             ( (cur_yr == end_yr) && (cur_per <= end_per)) ) {
        temp.insert(MISSING);
        if((cur_per%per_len) == 0){
          cur_per = 1;
          cur_yr++;
        }
        else
          cur_per++;
      }
    } 

  }
  else {
    if(start_yr == year[0])
      num_missing = (int)period[0] - start_per;
    else{
      if(period[0] <= start_per)
        num_missing = ((int)year[0] - start_yr - 1)*per_len + ((int)period[0] - start_per + per_len);
      else
        num_missing = ((int)year[0] - start_yr)*per_len + ((int)period[0] - start_per);
    }
    printf("num_missing=%d\n",num_missing);
    for(j = 0; j < num_missing; j++)
      temp.insert(MISSING);

    i = 0;
    while( ( ( year[i] < end_yr ) ||
                   ( year[i] == end_yr && period[i] <= end_per ) ) &&
                 ( i < List.get_size() ) ) {
        temp.insert(Array[i]);
 	printf("i=%d  cur_date: %d/%d  end_date: %d/%d\n",i,(int)period[i],(int)year[i],end_per,end_yr);
        i++;
    }
    if(i == List.get_size()){
      cur_yr = (int)year[i-1];
      cur_per = (int)period[i-1];
      if((cur_per%per_len) == 0){
        cur_per = 1;
        cur_yr++;
      }
      else
        cur_per++;   
      while( (cur_yr < end_yr) ||
             ( (cur_yr == end_yr) && (cur_per <= end_per)) ) {
        temp.insert(MISSING);
        printf("here i=%d  cur_date: %d/%d  end_date: %d/%d\n",i,cur_per,cur_yr,end_per,end_yr); 
        if((cur_per%per_len) == 0){
          cur_per = 1;
          cur_yr++;
        }
        else
          cur_per++;
      }
    }
  }
  delete [] Array;
  delete [] year;
  delete [] period;
  size = temp.get_size();
  return temp.returnArray();
}
inline int pdsi::is_int(char *string,int length) {
  int err=1; 
  for(int i=0; i<length; i++) 
    if(!isdigit(string[i])) err=0;
  return err; 
}
inline int pdsi::is_flt(char *string,int length) {
  int err=1;
  for(int i=0; i<length; i++)
    if(!isdigit(string[i]) && string[i]!='.') err=0;
  return err;
}
llist::llist() {
  head = new node;       // Creates the sentinel with head pointing to it
  head->next = head;     // Sets the sentinels front pointer to itself
  head->previous = head; // Sets the sentinels rear pointer to itself
  size = 0;
}
llist::~llist() {
  node *mover;                // The temporary pointer to perform the work
  
  mover = head->next;         // mover starts at the node after the sentinel
  while (mover != head) {     // This loop occurs until mover has come complete
                              // circle and is pointing at the sentinel
    mover = mover->next;      // mover becomes the next node
    delete mover->previous;   // The previous node is then deleted
  }
  delete mover;               // Finally the sentinel is deleted
}
void llist::insert(number x) {
  node *inserter;      // A new pointer to the node to be added
  inserter = new node; // This creates the node and points inserter to it
  inserter->key = x;   // The value in inserter is set to x
  inserter->next = head->next;        // The next field of the new node is
                                      // set to the node after the sentinel
  inserter->previous = head;          // The previous field is set to the 
                                      // sentinel
  inserter->next->previous = inserter;// The previous field of the node after
                                      // inserter is set to inserter
  inserter->previous->next = inserter;// The sentinels next field is set to 
                                      // inserter
  size++;                             // update size
}
int llist::get_size() {
  return size;
}
number* llist::returnArray() {
  node* cur;
  int i;
  number* A = new number[size];
  if(A != NULL){
	cur = head->previous;
	i = 0;
	while(cur != head){
	  A[i] = cur->key;
	  i++;
	  cur=cur->previous;
	}
  }
  return A;

}
number llist::head_remove() {
  if(is_empty()) {
    return MISSING;
  }
  node *remover;
  number x;
  remover=head->next;
  x=remover->key;
  // First the previous field of the next node is set to head
  remover->next->previous = head;
  // Then the next field of head is set to the node after remover
  head->next = remover->next;
  size--;          //update size;
  delete remover;  // Finally the node can be deleted and function ends
  return x;  // The key is returned
}
number llist::tail_remove() {
  if(is_empty()) {
    return MISSING;
  }
  node *remover;
  number x;
  
  remover=head->previous;
  x=remover->key;
  //First the next field of the previous node is set to head
  remover->previous->next = head;
  //Then the previous field of head is set to the node before remover
  head->previous = remover->previous;
  size--;            //update size
  delete remover;    // Finally the node can be deleted
  return x; // The key is returned
}
node *llist::set_node(node *set, number x) {
  int error=1;
  node *comparer;
  
  if (set==NULL)
    set=head->next;
  comparer = head->next;
  while(comparer != head) {
    if(comparer == set) {
      error=0;
      break;
    }
    comparer = comparer->next;
  }
  
  if(error==1) {
    return NULL;
  }
  else {
    if(set->key != MISSING){
      set->key = x;
      return set->next;
    }
    else {
      return set_node(set->next,x);
    }
  }
}
int llist::is_empty() {
  if(head->next==head)
    return 1;
  else
    return 0;
}
void copy(llist &L1,const llist &L2) {
  while (!L1.is_empty()) 
    L1.tail_remove();
  node *comparer;
  comparer = L2.head->previous;
  while (comparer!=L2.head) {
    L1.insert(comparer->key);
    comparer = comparer->previous;
  }
}
number llist::sumlist(){
  number sum = 0;
  node* cur;
  cur = head->previous;
  while(cur != head){
    sum += cur->key;
    cur = cur->previous;
  }

  return sum;
}
void llist::sumlist(number &prev_sum, int sign){
  number sum = 0;
  node* cur;
  cur = head->previous;
  while(cur != head){
    sum += cur->key;
    cur = cur->previous;
  }
  if(sign*sum > sign * prev_sum)
    prev_sum = sum;
  return;
}
number llist::maxlist(){
  number max = 0;
  node * cur;
  cur = head->previous;
  while(cur != head){
    if(cur->key > max)
      max = cur->key;
    cur = cur->previous;
  }
  return max;
}
number llist::minlist(){
  number min = 0;
  node * cur;
  cur = head->previous;
  while(cur != head){
    if(cur->key < min)
      min = cur->key;
    cur = cur->previous;
  }

  return min;
}
number llist::kthLargest(int k) {
  if(k < 1 || k > size)
    return MISSING;
  else if(k == 1)
    return minlist();
  else if(k == size) 
    return maxlist();
  else{
    //place the list in an array for
    //easier selection
    number *A;
    int i; 
    number return_value;
    node* cur = head->previous;
    A = new number[size];
    if(A != NULL){
      for(i = 0; i < size; i++){
		if(cur->key != MISSING)
		  A[i] = cur->key;
		cur = cur->previous;
      }
      select(A,0,size-1,k);
      
      return_value = A[k-1];
      delete []A;
    }
    else {
      long_select(k);
    }
    return return_value;
  }
}
number llist::quartile(int q) {
  //q0 is the minimum
  if(q == 0)
    return minlist();
  //q4 is the maximum
  else if(q == 4)
    return maxlist();
  else if(q == 2) {
    if(size%2 == 0){
      double t1 = kthLargest(size/2);
      double t2 = kthLargest(size/2 + 1);
      return (t1+t2)/2;
    }
    else
      return kthLargest(1 + (size-1)/2);
  }
  
  //q1 is the first quartile, q3 is the third quartile
  else if(q == 1 || q == 3){
    //if (size+1) is not divisble by 4, there is no exact quartiles
    //so take the weighted average of the two closest numbers
    int k;
    if((k = ((size-1)%4)) != 0){
      int bottom = (int)floor(q*(size-1)/4);
      double t1 = (4-k) * kthLargest(bottom+1);
      double t2 = (k) * kthLargest(bottom+2);
      return (t1+t2)/4;
    }
    else
      return kthLargest(1 + q*(size-1)/4);
  }
  else
    return MISSING;
}
number llist::safe_percentile(double percentage) {
  llist temp;
  node* cur = head->next; 
  while(cur != head){ 
    if(cur->key != MISSING) 
      temp.insert(cur->key);
    cur = cur->next; 
  } 
  return temp.percentile(percentage);
}
number llist::percentile(double percentage) {
  int k;
  if(percentage > 1)
    percentage = percentage / 100;
  k = (int)(percentage * size);
  return kthLargest(k);
}
number llist::long_select(int k) {
  printf("Low Memory.\n");
  return MISSING;
  
}
int numEntries(FILE *in) {
  int i = 0;
  float t;
  while(fscanf(in,"%f",&t)!=EOF)
    i++;
  return i;
}
int partition(number a[], int left, int right) {
  number val = a[right];
  int i = left - 1;
  for(int j = left; j < right; j++){
	if(a[j] <= val){
	  i++;
	  exch(a[i],a[j]);
	}
  }
  exch(a[i+1],a[right]);
  return i+1;
}
void select(number a[], int l, int r, int k) {
  int i;
  if (r <= l) 
	return;

  i = partition(a, l, r);
  if (i > k-1) 
	select(a, l, i-1, k);
  else if (i < k-1) 
	select(a, i+1, r, k);
  else 
	return;
}
void exch(number &x, number &y) {
  number temp;
  temp = x; 
  x = y;
  y = temp;
}
int dir_exists(char *dir) {
  FILE* test;
  char test_file[128];
  if(dir[strlen(dir)-1] != '/')
    strcat(dir,"/");
  strcpy(test_file,dir);
  strcat(test_file, "test.file");
  test = fopen(test_file, "w");
  if(test == NULL){
    return -1;
  }
  fclose(test);
  remove(test_file);
  return 1;
}
int create_dir(char *path) {
  unsigned int i = 0;
  char my_path[128];
  char dir[128];
  int return_value = 1;
  strcpy(my_path,path);
  if(path[strlen(path)-1] != '/')
    strcat(my_path,"/");

  while(i < strlen(my_path)){
    dir[i] = my_path[i];
    if(my_path[i] == '/'){
      dir[i+1] = '\0';
      //check to see if this directory exists
      if(dir_exists(dir) == -1){
		return_value = 0;
		if(_mkdir(dir)!=0)
		  return -1;
      }
    }
    i++;
  }
  return return_value;
  //---------------------------------------------------------------------------
}
