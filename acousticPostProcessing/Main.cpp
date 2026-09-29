#include "FileIOTools.h"
#include "FFT.h"
#include "json.hpp"

int main(int argc, char** argv)
{
  assert(argc == 2);
  std::string jsonPath = argv[1];
  std::ifstream jsonFile(jsonPath);
  nlohmann::json data = nlohmann::json::parse(jsonFile);
  double minFreq = *(data.find("minimumFrequency"));
  double maxFreq = *(data.find("maximumFrequency"));
  std::vector<std::string> filePathVec = *(data.find("filePath"));
  int windowFuncType = *(data.find("windowFunctionType"));
  std::string band = *(data.find("band"));
  std::string weighting = *(data.find("weighting"));
  double startTime = *(data.find("startTime"));
  double endTime = *(data.find("endTime"));
  std::string writePath = *(data.find("writePath"));
  int nSmooth = 5;
  if (data.contains("smoothFactor"))
  {
    nSmooth = *(data.find("smoothFactor"));
  }

  //* 检查输入的内容
  if (!(windowFuncType >= 1 && windowFuncType <= 6))
  {
    std::cout << "Incorrect window function type." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(band == "narrowBand" || band == "oneThirdOctaveBand"))
  {
    std::cout << "Incorrect band." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(weighting == "A-weight" || weighting == "none"))
  {
    std::cout << "Incorrect weighting." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (filePathVec.empty())
  {
    std::cout << "No input file." << std::endl;
    exit(EXIT_SUCCESS);
  }
  for (const auto& filePath : filePathVec)
  {
    if (!std::ifstream(filePath.c_str()).good())
    {
      std::cout << "File \"" << filePath << "\" not found." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //* 从每个文件第一行获取测点坐标
  std::vector<std::array<double, 3>> pointPosVec;
  for (const auto& filePath : filePathVec)
  {
    pointPosVec.push_back(read_point_position(filePath));
  }

  //* 读取时间轴、时域的近场测点压力/远场测点声压
  std::vector<std::vector<double>> timeVec;
  std::vector<std::vector<double>> preTimeDomain;
  read_from_file(timeVec, preTimeDomain, filePathVec, startTime, endTime);

  //* 快速傅里叶变换：得到频率轴(>=minFreq, <=maxFreq)、频域的声压
  std::vector<std::vector<double>> freqVec;
  std::vector<std::vector<std::complex<double>>> preFreqDomain;
  cal_fft(freqVec, preFreqDomain, timeVec, preTimeDomain, windowFuncType);

  //* 考虑倍频程(band)、计权方式(weighting)，频域的声压计算频域的声压级
  std::vector<std::vector<double>> freqWeightedBandVec;
  std::vector<std::vector<double>> splWeightedBandVec;
  std::vector<std::vector<double>> freqSmoothedVec;
  std::vector<std::vector<double>> splWeightedSmoothedVec;
  std::vector<double> oaSplWeightedSmoothedVec = pressure_to_spl_smoothed(freqWeightedBandVec, splWeightedBandVec, freqSmoothedVec, splWeightedSmoothedVec, freqVec, preFreqDomain, minFreq, maxFreq, band, weighting, nSmooth);

  //* 由计权的声压级计算功率谱密度级
  std::vector<std::vector<double>> psdWeightedSmoothedVec;
  std::vector<std::vector<double>> psdWeightedSmoothedLevelVec;
  spl_to_psd_smoothed(psdWeightedSmoothedVec, psdWeightedSmoothedLevelVec, freqSmoothedVec, splWeightedSmoothedVec);

  // //* 考虑频率范围准备输出
  // for (size_t iObs = 0; iObs < freqVec.size(); ++iObs)
  // {
  //   apply_freq_limits(freqWeightedBandVec[iObs], splWeightedBandVec[iObs], freqWeightedBandVec[iObs], splWeightedBandVec[iObs], minFreq, maxFreq);
  //   apply_freq_limits(freqSmoothedVec[iObs], splWeightedSmoothedVec[iObs], freqWeightedBandVec[iObs], splWeightedBandVec[iObs], minFreq, maxFreq);
  //   apply_freq_limits(freqSmoothedVec[iObs], psdWeightedSmoothedVec[iObs], freqSmoothedVec[iObs], psdWeightedSmoothedVec[iObs], minFreq, maxFreq);
  //   apply_freq_limits(freqSmoothedVec[iObs], psdWeightedSmoothedLevelVec[iObs], freqSmoothedVec[iObs], psdWeightedSmoothedLevelVec[iObs], minFreq, maxFreq);
  // }

  //* 将各测点的声压级（考虑计权、倍频程）写到文件，一个测点对应一个文件
  write_SPL_or_PSD_to_file(freqWeightedBandVec, splWeightedBandVec, pointPosVec, writePath, weighting, band, "SPL");
  //* 将各测点的声压级（考虑计权）写到文件，一个测点对应一个文件
  write_SPL_or_PSD_to_file(freqSmoothedVec, splWeightedSmoothedVec, pointPosVec, writePath, weighting, "narrowBand", "SPL");
  //* 将各测点的功率谱密度级（考虑计权）写到文件，一个测点对应一个文件
  write_SPL_or_PSD_to_file(freqSmoothedVec, psdWeightedSmoothedLevelVec, pointPosVec, writePath, weighting, band, "PSDL");

  //* 打印各测点的总声压级
  for (size_t iObs = 0; iObs < freqVec.size(); ++iObs)
  {
    std::cout << "iObs " << iObs << " OASPL = " << oaSplWeightedSmoothedVec[iObs] << "dB(A)" << std::endl;
  }

  return 0;
}
