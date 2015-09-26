//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "path_helper.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>

namespace path {

bool make_dir(std::string const& folder) {
    return mkdir(folder.c_str(), 0755) == 0;
}

std::string make_temp_dir(std::string const& prefix,
                          std::string const& suffix) {
    std::string name = append_path(prefix, suffix + "_XXXXXX");
    char* actual;
    if ((actual = ::mkdtemp(strcpy(new char[name.length() + 1], name.c_str())))
            == NULL)
        throw std::runtime_error("Cannot create temporary dir " + name);

    std::string result(actual);
    if (result == name)
        throw std::runtime_error("Cannot create temporary dir " + name);

    delete[] actual;

    return result;
}

void remove_dir(std::string const& folder) {
    DIR *dp;
    if ((dp = opendir(folder.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + folder);

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL) {
        std::string full_path = folder + "/" + dirp->d_name;

        if (dirp->d_type == DT_DIR) {
            if (std::string(".") != dirp->d_name
                    && std::string("..") != dirp->d_name) {
                remove_dir(full_path);
            }
        } else
            remove(full_path.c_str());
    }

    closedir(dp);
    remove(folder.c_str());
}

bool is_regular_file(std::string const& path) {
    struct stat st;
    return (stat(path.c_str(), &st) == 0) && (S_ISREG(st.st_mode));
}

std::string append_path(std::string const& prefix, std::string const& suffix) {
    std::string delimiter = "";

    if (!boost::ends_with(prefix, "/") && !boost::starts_with(suffix, "/")
            && !prefix.empty()) {
        delimiter = "/";
    }

    return prefix + delimiter + suffix;
}

std::string current_dir() {
    char* cwd = getcwd(NULL, 0);
    std::string result = cwd;

    free(cwd);
    return result;
}

void make_full_path(std::string& path) {
    if (!boost::starts_with(path, "/"))  // relative path
        path = append_path(current_dir(), path);
}

std::string filename(std::string const& path) {
    size_t pos = path.find_last_of('/');
    return pos != std::string::npos ? path.substr(pos + 1) : path;
}

std::string basename(std::string const& path) {
    size_t slash = path.find_last_of('/');
    size_t after_slash = slash == std::string::npos ? 0 : slash + 1;

    size_t dot = path.find_last_of('.');
    if (dot < after_slash)
        dot = std::string::npos;

    return path.substr(after_slash, dot - after_slash);
}

std::string extension(std::string const& path) {
    size_t slash = path.find_last_of('/');
    size_t after_slash = slash == std::string::npos ? 0 : slash + 1;
    size_t dot = path.find_last_of('.');

    if (dot < after_slash || dot == std::string::npos || dot + 1 == path.size())
        return std::string();

    return path.substr(dot);
}

std::string parent_path(std::string const& path) {
    std::string cpath(path);

    make_full_path(cpath);
    size_t slash_pos = cpath.find_last_of('/');

    return (slash_pos == 0 ? std::string("/") : cpath.substr(0, slash_pos));
}

bool check_existence(std::string const& path) {
    struct stat st_buf;
    return stat(path.c_str(), &st_buf) == 0
            && (S_ISREG(st_buf.st_mode) || S_ISDIR(st_buf.st_mode));  // exists and (file or dir)
}

void remove_if_exists(std::string const& path) {
    if (check_existence(path)) {
        if (is_regular_file(path)) // file
            remove(path.c_str());
        else // dir
            remove_dir(path);
    }
}

// doesn't support symlinks
std::string resolve(std::string const& path) {
    typedef boost::char_delimiters_separator<char> separator_t;
    typedef boost::tokenizer<separator_t> tokenizer_t;

    tokenizer_t tok(path, separator_t(false, "", "/"));

    std::string result = "/";
    for (auto it = tok.begin(); it != tok.end(); ++it) {
        if (*it == "..")
            result = parent_path(result);

        else if (*it == ".")
            ;  // Ignore

        else
            // Just cat other path entries
            result = append_path(result, *it);
    }

    return result;
}

std::string make_relative_path(std::string p, std::string base) {
    p = resolve(p);
    base = resolve(base);

    std::string pp = parent_path(p);

    typedef boost::char_delimiters_separator<char> separator_t;
    typedef boost::tokenizer<separator_t> tokenizer_t;

    tokenizer_t pp_tok(pp, separator_t(false, "", "/"));
    tokenizer_t base_tok(base, separator_t(false, "", "/"));

    auto i = pp_tok.begin();
    auto j = base_tok.begin();

    while (i != pp_tok.end() && j != base_tok.end() && *i == *j) {
        ++i;
        ++j;
    }

    std::string result;
    for (; j != base_tok.end(); ++j)
        result = append_path("..", result);

    for (; i != pp_tok.end(); ++i)
        result = append_path(result, *i);

    return append_path(result, filename(p));
}

typedef std::vector<std::string> files_t;

}
