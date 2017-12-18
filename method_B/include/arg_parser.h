#pragma once

#include <cassert>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>

/// Parse command-line arguments
/**
 * A simple command-line parser.
 *
 * Supports named arguments and switches as well as unnamed data arguments
 *
 * Example: ./foo -v -o outfolder in1.xml in2.xml
 */
class arg_parser {
public:
    /// Parse command-line arguments
    arg_parser(int argc, char **argv) : named_args(), data_args() {
        int pos = 1;
        while (pos < argc) {
            if (argv[pos][0] == '-') {
                // Read argument and advance
                std::string arg(&argv[pos++][1]);
                // check whether argument has a value (i.e. not just a flag)
                if (pos < argc && argv[pos][0] != '-') {
                    named_args[arg] = std::string{argv[pos++]};  // assign and advance to next
                } else {
                    named_args[arg] = "";
                }
            } else {
                data_args.emplace_back(argv[pos++]);
            }
        }
    }

    /// Get a named argument's value
    /// \param key the argument name
    /// \param default_value the value to return if the argument wasn't set
    template <typename T>
    T get(const std::string &key, const T default_value = T()) const {
        T retval;
        auto it = named_args.find(key);
        if (it != named_args.end()) {
            std::istringstream s(it->second);
            s >> retval;
        } else {
            // do this in the else case, otherwise empty string arguments
            // would return the default value instead of ""
            retval = default_value;
        }
        return retval;
    }

    /// check whether an argument was set
    bool is_set(const std::string &arg) const {
        return named_args.find(arg) != named_args.end();
    }

    /// the number of unnamed data arguments
    size_t num_data_args() const {
        return data_args.size();
    }

    /// get a data argument by its index (among the data arguments)
    std::string data_arg(const size_t index) const {
        assert(index < num_data_args());
        return data_args[index];
    }

protected:
    std::unordered_map<std::string, std::string> named_args;
    std::vector<std::string> data_args;
};
