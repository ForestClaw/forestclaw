/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <fclaw_filesystem.h>
#include <fclaw_base.h>
#if __APPLE__ || __linux__
#include <unistd.h>
#else
#include <filesystem>
#include <string>
#include <cstring>
#endif

char* fclaw_cwd()
{
#if __APPLE__ || __linux__
    char* c_current_path = FCLAW_ALLOC(char,PATH_MAX+1);
    char* error = getcwd(c_current_path,PATH_MAX+1);
    FCLAW_ASSERT(error != NULL);
    return c_current_path;
#else
    std::string current_path = std::filesystem::current_path();
    char* c_current_path = FCLAW_ALLOC(char, current_path.length()+1);
    strcpy(c_current_path,current_path.c_str());
    return c_current_path;
#endif
}

void fclaw_cd(const char* dir)
{
#if __APPLE__ || __linux__
    int error = chdir(dir);
    FCLAW_ASSERT(error == 0);
#else
    std::filesystem::path output_dir(dir);
	std::filesystem::path output_path;
    if(output_dir.is_absolute())
    {
	    output_path = output_dir;
    }
    else
    {
        output_path = std::filesystem::current_path();
        output_path /= output_dir;
    }

	std::filesystem::create_directory(output_path);

	std::filesystem::current_path(output_path);
#endif
}

