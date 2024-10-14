#ifndef PATH_NAMES_DATA_CODE
#define PATH_NAMES_DATA_CODE 1

// const char* 
// #define DATA_PATH_NAME "/diskb/nBLM/data/production"
#define DATA_PATH_NAME "/eos/home-t/tpapaeva/picosec/data/2022_July_h4"
// const char*
#define WORK_DIR_NAME "/eos/home-t/tpapaeva/picosec/data/"
#define CODEDIR "/eos/home-t/tpapaeva/picosec/code/"

#define RUN_TYPE "PICOSECraw"

#define RUNMIN 1
#define RUNMAX 1500

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KRESET "\x1B[0m"


const std::string RED("\033[0;31m"); 
const std::string RED_U_GBKG("\033[4;31;42m"); 
const std::string GREEN("\033[0;32m"); 
const std::string YELLOW("\033[0;33m"); 
const std::string BLUE("\033[0;34m"); 
const std::string MAGENTA("\033[0;35m"); 
const std::string CYAN("\033[0;36m"); 
const std::string INVERSE_ON("\033[7m"); 
const std::string INVERSE_OFF("\033[27m"); 
const std::string RESET_COLOR("\033[0m");	
const std::string endlr("\n\033[0m");	


#endif
