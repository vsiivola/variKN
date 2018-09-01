#include "io.hh"

namespace io {
bool Stream::verbose = false;

Stream::Stream() : file(NULL), is_pipe(false), close_allowed(true) {}

Stream::Stream(std::string file_name, std::string mode, SeekType seek_type,
               bool allow_close)
    : file(NULL), is_pipe(false), close_allowed(true) {
  open(file_name, mode, seek_type, allow_close);
}

Stream::~Stream() { close(); }

void Stream::open(std::string file_name, std::string mode, SeekType seek_type,
                  bool allow_close) {
  close();

  if (file_name.length() == 0) {
    fprintf(stderr, "can not open empty filename\n");
    throw OpenError();
  }

  close_allowed = allow_close;

  if (verbose) {
    if (mode == "r")
      fprintf(stderr, "Opened read: ");
    else
      fprintf(stderr, "Opened write: ");
    fprintf(stderr, "%s\n", file_name.c_str());
  }

  // Standard input or output
  if (file_name == "-") {
    if (seek_type != LINEAR) {
      fprintf(
          stderr,
          "Stream::open(): "
          "Program requested reopenable or seekable input, but got %s. Exit\n",
          file_name.c_str());
      throw OpenError();
    }
    if (mode == "r") {
      file = stdin;
    } else if (mode == "w")
      file = stdout;
    else {
      fprintf(stderr,
              "Stream::open(): "
              "invalid mode %s for standard input or output\n",
              mode.c_str());
      throw OpenError();
    }
    close_allowed = false;
  }

  // A pipe for reading
  else if (file_name[file_name.length() - 1] == '|') {
    if (seek_type == SEEKABLE) {
      fprintf(stderr,
              "Stream::open(): "
              "Program requested seekable input, but got %s. Exit\n",
              file_name.c_str());
      throw OpenError();
    }
    /*if (mode != "r") {
      fprintf(stderr, "invalid mode %s for input pipe %s\n",
              mode.c_str(), file_name.c_str());
      exit(1);
      }*/
    file_name.erase(file_name.length() - 1, 1);
#ifdef _WIN32
    throw OpenError();
#else
    file = popen(file_name.c_str(), mode.c_str());
#endif
    if (file == NULL) {
      fprintf(stderr, "could not open pipe %s\n", file_name.c_str());
      throw OpenError();
    }
    is_pipe = true;
  }

  // A pipe for writing
  else if (file_name[0] == '|') {
    if (seek_type == SEEKABLE) {
      fprintf(stderr,
              "Stream::open(): "
              "Program requested seekable input, but got %s. Exit\n",
              file_name.c_str());
      throw OpenError();
    }
    /*if (mode != "w") {
      fprintf(stderr, "invalid mode %s for output pipe %s\n",
              mode.c_str(), file_name.c_str());
      exit(1);
      }*/
    file_name.erase(0, 1);
#ifdef _WIN32
    throw OpenError();
#else
    file = popen(file_name.c_str(), mode.c_str());
#endif
    if (file == NULL) {
      fprintf(stderr, "could not open pipe %s\n", file_name.c_str());
      throw OpenError();
    }
    is_pipe = true;
  }

  // A gzipped file
  else if ((file_name.length() >= 3 &&
            file_name.substr(file_name.length() - 3) == ".gz") ||
           (file_name.length() >= 4 &&
            file_name.substr(file_name.length() - 4) == ".bz2"))

  {
    if (seek_type == SEEKABLE) {
      fprintf(stderr,
              "Stream::open(): "
              "Program requested seekable input, but got %s. Exit\n",
              file_name.c_str());
      throw OpenError();
    }
    std::string pipe;
    if (mode == "r") {
      if (file_name.substr(file_name.length() - 3) == ".gz")
        pipe = "gzip -dc '";
      else
        pipe = "bzip2 -dc '";
      pipe += file_name + "'";
    } else if (mode == "w") {
      if (file_name.substr(file_name.length() - 3) == ".gz")
        pipe = "gzip > '";
      else
        pipe = "bzip2 --best > '";
      pipe += file_name + "'";
    } else {
      fprintf(stderr, "invalid mode %s for gzipped file\n", mode.c_str());
    }
#ifdef _WIN32
    throw OpenError();
#else
    file = popen(pipe.c_str(), mode.c_str());
#endif
    if (file == NULL) {
      fprintf(stderr, "could not open mode %s for pipe: %s\n", mode.c_str(),
              pipe.c_str());
      throw OpenError();
    }
    is_pipe = true;
  }

  // An ordinary file
  else {
    file = fopen(file_name.c_str(), mode.c_str());
    if (file == NULL) {
      fprintf(stderr, "could not open mode %s for file %s\n", mode.c_str(),
              file_name.c_str());
      throw OpenError();
    }
    is_pipe = false;
  }
}

void Stream::close() {
  if (file == NULL || !close_allowed)
    return;
  if (is_pipe)
#ifdef _WIN32
    throw OpenError();
#else
    pclose(file);
#endif
  else
    fclose(file);
  is_pipe = false;
  file = NULL;
}

} // namespace io
