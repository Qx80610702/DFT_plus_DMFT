#include "utilities.h"

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <stdio.h>

int rm_dir(const char *dir)
{
  DIR *dir_tmp;
  struct dirent  *dirp;
  struct stat  buf;
  char* p = getcwd( NULL, 0 );

  if(access(dir,0) != 0) return 0;  //Directory does not exist

  if((dir_tmp = opendir(dir) ) == NULL ) return -1;

  chdir(dir);

  while( dirp = readdir(dir_tmp) )
  {
    if((strcmp( dirp->d_name, "." ) == 0) || (strcmp( dirp->d_name, ".." ) == 0))
      continue;

    if(stat( dirp->d_name, &buf ) == -1) return -1;

    if( S_ISDIR(buf.st_mode) ){
      rm_dir( dirp->d_name );
      /*if(rmdir(dirp->d_name)==-1)
       *  error_quit("rmdir");
       * printf("rm %s Successed . . .\n",dirp->d_name);*/
      continue;
    }

    if(remove(dirp->d_name) == -1) return -1;
  }
  closedir(dir_tmp);
  chdir(p);

  if (rmdir(dir) == -1) return -1;

  return 0;
}

int mv_dir(const char *dir_old, const char *dir_new)
{
  return rename(dir_old, dir_new);
}


int mk_dir(const char *dir)
{
  return mkdir(dir,S_IRWXU);
}

int cp_dir(const char *dir_old, const char *dir_new)
{
  
}
