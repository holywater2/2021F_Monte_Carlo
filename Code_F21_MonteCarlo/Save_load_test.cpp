#include <sys/types.h>
#include <direct.h>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

bool file_exists(string filename) {
    struct stat fileinfo;

    return !stat(filename.c_str(), &fileinfo);
}

void save_to_txt(string _data, string _filename){
    int idx = 0;
    string folder_name = make_next_directory("result");
    
    string filePath = folder_name + string("/test.txt");

	// write File
	ofstream writeFile(filePath.data());
	if( writeFile.is_open() ){
		writeFile << "Hello World!\n";
		writeFile << "This is C++ File Contents.\n";
		writeFile.close();
	}

	// read File
	ifstream openFile(filePath.data());
	if( openFile.is_open() ){
		string line;
		while(getline(openFile, line)){
			cout << line << endl;
		}
		openFile.close();
	}
}

// Make New directory for save data
// ex) input: result
//      if there is result0, then make "result1" folder.
//     return result1
string make_next_directory(string folder_name){
    int idx = 0;
    while(true){
        if( !file_exists(folder_name+char(idx+'0')) ){
            break;
        }
        idx++;
    }
    folder_name += char(idx+'0');
    mkdir(folder_name.c_str());
    return folder_name;
}