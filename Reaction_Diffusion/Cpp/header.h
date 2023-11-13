#pragma once

#include <iostream>
#include <fstream>
#include <boost/chrono.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <stdio.h>

#include "ReactionDiffusion.h"

// input-handler declaration
void input_handler(int& argc, char** argv, ReactionDiffusion& System);
