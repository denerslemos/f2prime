#include "call_libraries.h"  // call libraries from ROOT and C++

// print UIC Jets welcome message
/*
Arguments
printout: true print UIC Jet logo, false do not print out
*/
void printwelcome(bool printout){
	if(!printout) return;
	cout << endl;
	cout << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"         << endl;
	cout << "+                                                              +"         << endl;
	cout << "+    	        Welcome to K0Star skim code for Data     	+"         << endl;
	cout << "+              From Dener Lemos: ddesouzal@bnl.gov             +"         << endl;
	cout << "+                                                              +"         << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"         << endl;
}

// print start day and time message
void print_start(){
	cout << endl;
    time_t init = time(0);
    char* init_time = ctime(&init); // convert now to string form
    cout << "Starting at : " << init_time << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
	cout << endl;
}

// print stop day and time message
void print_stop(){
    time_t end = time(0);
    char* end_time = ctime(&end); // convert now to string form
   	cout << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   	cout << endl;
	cout << "Stopping at : " << end_time << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << endl;
}
