
#ifndef _agb_model_h
#define _agb_model_h

#include <iostream>
#include <vector>
#include <string>

#include "../config/config.h"
#include "../spectral_grid/spectral_grid.h"

namespace agb {


class AGBStarModel{
  public:
    AGBStarModel(const std::string folder);
    ~AGBStarModel() {}

    void calcModel();

    ModelConfig config;
    SpectralGrid spectral_grid;
  protected:
};



}

#endif