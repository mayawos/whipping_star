#include <getopt.h>

#include "SBNcovariance.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{

  std::string xml = "";
  int iarg = 0;
  opterr = 1;
  int index;

  /*************************************************************
   *************************************************************
   *		Command Line Argument Reading
   ************************************************************
   ************************************************************/
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"test",		required_argument,	0, 't'},
      {0,			no_argument, 		0,  0},
    };

  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 't':
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }

  /*************************************************************
   *************************************************************
   *			Main Program Flow
   ************************************************************
   ************************************************************/
  time_t start_time = time(0);
  std::cout<<"Begining Covariance Calculation: "<<std::endl;

  //a tag to identify outputs
  std::string tag = "DL";

  //Create a SBNmultisim object initilizing with the inputted xml
  SBNcovariance lee_covar(xml);

  //Form the covariance matrix from loaded weights
  lee_covar.FormCovarianceMatrix(tag);

  //and make some plots of the resulting things
  lee_covar.PrintMatricies(tag);

  std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
  return 0;

}
