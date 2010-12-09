#ifndef OUTPUT_H
#define OUTPUT_H
#include <boost/filesystem.hpp>
#include <string>
namespace fs = boost::filesystem;
using namespace std;

class t_output
{
    string output_dir;
    public:
    t_output(string _output_dir) : output_dir(_output_dir)
    {
	fs::path output_path(output_dir);
	if( fs::exists(output_path) )
	{
	    cout << output_path.leaf() <<" already exists, use it anyway? y/n\n";
	    string response;
	    cin >> response;
	    if(response!="y") exit(1);
	    if(!fs::is_directory(output_path))
	    {
		cerr << output_path.leaf() << " is not a directory\n";
		exit(1);
	    }
	}
	else
	{
	    cout << output_path << endl;
	    fs::create_directories(output_path);
	}
    }
    void backup(string oldname, string newname)
    {
        // make backup of configuration
        string config_backup = output_dir + "/" + newname;
        cout << "making backup of " + oldname + " to " + config_backup << endl;
        if(fs::exists(config_backup) && !fs::is_directory(config_backup))
            fs::remove(config_backup);
        fs::copy_file(oldname,config_backup);

    }

};


#endif
