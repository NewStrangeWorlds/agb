

#include <iostream>
#include <omp.h>
#include <csignal>
#include <cstdlib>

#include "../agb_model/agb_model.h"


int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "wrong number of command line parameters!\n";

    return 1;
  }

  std::string model_folder = argv[1];


  agb::AGBStarModel agb_wind(model_folder);

  agb_wind.calcModel();

  return 0;
}