#include <array>
#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <filesystem>

/** @brief 从结果文件第一行读测点坐标
 * @param[in] filePath 文件路径
 * @return std::array<double, 3> 测点坐标
 */
std::array<double, 3> read_point_position(const std::string& filePath)
{
  std::ifstream file(filePath, std::ios::in);

  //* first line
  std::string firstLine;
  std::getline(file, firstLine);

  // remove parentheses
  firstLine.erase(0, 12);                    // remove the opening parenthesis
  firstLine.erase(firstLine.size() - 1, 1); // remove the closing parenthesis

  // create a stringstream to process the firstLine
  std::stringstream ssPt(firstLine);
  std::vector<double> pt;
  double number;
  while (std::getline(ssPt, firstLine, ','))
  {
    std::stringstream num(firstLine);
    num >> number;
    pt.push_back(number);
  }
  file.close();
  return std::array<double, 3>{pt[0], pt[1], pt[2]};
}

/** @brief 读文件第二行各列标题
 * @param[in] filePath 文件路径
 * @return std::vector<std::string>第二行标题
 */
std::vector<std::string> read_titleVec(const std::string& filePath)
{
  std::ifstream file(filePath, std::ios::in);
  //* second line
  std::string secondLine;
  std::getline(file, secondLine);
  std::getline(file, secondLine);
  std::stringstream ssTitle(secondLine);
  std::vector<std::string> titleVec;
  std::string titleStr;
  while (ssTitle >> titleStr)
  {
    titleVec.push_back(titleStr);
  }
  file.close();
  return titleVec;
}

/** @brief 从文件获取指定标题对应的列内容
 * @param[in] filePath 文件路径
 * @param[in] skipHeader 跳过文件头行数
 * @param[in] columnTitle 需获取列内容的标题
 * @param[in] titleVec 所有列标题
 * @return 标题对应的列数据
 */
template<typename T>
std::vector<T> read_column(
  const std::string& filePath,
  const int& skipHeader,
  const std::string& columnTitle,
  const std::vector<std::string>& titleVec)
{
  std::ifstream file(filePath, std::ios::in);
  std::string line;
  for (int countLine = 0; countLine < skipHeader; ++countLine)
    std::getline(file, line);

  int targetColumnIdx = 0;
  for (const auto& title : titleVec)
  {
    if (title == columnTitle)
      break;
    else
      targetColumnIdx++;
  }
  std::vector<T> dataVec;
  T dataTmp;

  while (std::getline(file, line))
  {
    std::stringstream ssData(line);
    int idx = 0;
    while (ssData >> dataTmp)
    {
      if (idx == targetColumnIdx)
        dataVec.push_back(dataTmp);
      idx++;
    }
  }

  return dataVec;
}

/** @brief 从所有测点文件读取时间轴、时域下的声压数据
 * @param[out] timeVec 不同测点的时间轴 obsNum signalLength
 * @param[out] preTimeDomain 时域下不同测点的声压 obsNum signalLength
 * @param[in] filePathVec 不同测点对应的文件路径
 * @param[in] startTime 开始时间
 * @param[in] endTime 结束时间
 */
void read_from_file(
  std::vector<std::vector<double>>& timeVec,
  std::vector<std::vector<double>>& preTimeDomain,
  const std::vector<std::string>& filePathVec,
  const double& startTime,
  const double& endTime)
{
  size_t obsNum = filePathVec.size();
  std::string title = "Pressure";
  timeVec.resize(obsNum);
  preTimeDomain.resize(obsNum);
  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    //* timeVec
    std::vector<std::string> titleVec = read_titleVec(filePathVec[iObs]);
    auto timeVecFile = read_column<double>(filePathVec[iObs], 2, "Time", titleVec);
    size_t signalLength = timeVecFile.size();
    size_t startTimeIdx = 0;
    size_t endTimeIdx = 0;
    double dt = (timeVecFile.back() - timeVecFile[0]) / double(signalLength - 1);
    startTimeIdx = size_t((startTime - timeVecFile[0]) / dt) + 1;
    endTimeIdx = size_t((endTime - timeVecFile[0]) / dt) + 1;
    if (startTime < timeVecFile[0])
    {
      std::cout << "测点编号" << iObs << "时域起始时间需大于" << timeVecFile[0] << "s. 已自动调整为" << timeVecFile[0] << "s." << std::endl;
      startTimeIdx = 0;
    }
    if (endTime > timeVecFile[signalLength - 1])
    {
      std::cout << "测点编号" << iObs << "时域终止时间需小于" << timeVecFile[signalLength - 1] << "s. 已自动调整为" << timeVecFile[signalLength - 1] << "s." << std::endl;
      endTimeIdx = signalLength - 1;
    }
    if (endTimeIdx >= signalLength)
      endTimeIdx = signalLength - 1;

    std::vector<double>::const_iterator startIt = timeVecFile.begin() + startTimeIdx;
    std::vector<double>::const_iterator endIt = timeVecFile.begin() + endTimeIdx + 1;
    timeVec[iObs].assign(startIt, endIt);
    // std::cout << "timeVec[0]: " << timeVec[0] << ", timeVec.back(): " << timeVec.back() << std::endl;

    //* pressureVec
    preTimeDomain[iObs] = read_column<double>(filePathVec[iObs], 2, title, titleVec);
    startIt = preTimeDomain[iObs].begin() + startTimeIdx;
    endIt = preTimeDomain[iObs].begin() + endTimeIdx + 1;
    preTimeDomain[iObs].assign(startIt, endIt);
  }
}

/** @brief 将测点对应的声压级or功率谱密度(级)输出到文件
 * @param[in] freqVec 频率轴
 * @param[in] splFreqDomain 频域声压级
 * @param[in] pointPosVec 测点坐标
 * @param[in] writePath 输出路径
 * @param[in] weighting 计权方式
 */
inline void write_SPL_or_PSD_to_file(
  const std::vector<std::vector<double>>& freqVec,
  const std::vector<std::vector<double>>& splOrPsdVec,
  const std::vector<std::array<double, 3>>& pointPosVec,
  const std::string& writePath,
  const std::string& weighting,
  const std::string& band,
  const std::string& splOrPsdOrPsdl)
{
  assert(splOrPsdVec.size() == pointPosVec.size());
  size_t obsNum = splOrPsdVec.size();
  std::string writePathCommon = writePath + (writePath.back() == '/' ? "" : "/");
  std::string aweight = (weighting == "A-weight" ? "A" : "");
  std::string unit;
  if (splOrPsdOrPsdl == "SPL")
  {
    if (band == "oneThirdOctaveBand")
    {
      unit = "SPL(1/3)(dB" + aweight + ")";
      writePathCommon += "sound_pressure_level_octave_pt";
    }
    else if (band == "narrowBand")
    {
      unit = "SPL(dB" + aweight + ")";
      writePathCommon += "sound_pressure_level_pt";
    }
  }
  else if (splOrPsdOrPsdl == "PSD")
  {
    unit = "PSD(dB" + aweight + "/Hz)";
    writePathCommon += "power_spectrum_density_pt";
  }
  else if (splOrPsdOrPsdl == "PSDL")
  {
    unit = "PSDL(dB" + aweight + "/Hz)";
    writePathCommon += "power_spectrum_density_level_pt";
  }
  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    std::stringstream ss;
    ss.setf(std::ios::right);
    ss << std::setw(4) << std::setfill('0') << std::to_string(iObs);
    std::string writePathFull = writePathCommon + ss.str() + ".txt";
    std::filesystem::create_directories(writePath);
    std::ofstream outFile(writePathFull, std::ios::out);
    outFile.setf(std::ios::left | std::ios::fixed);
    outFile << "(" << pointPosVec[iObs][0] << ", " << pointPosVec[iObs][1] << ", " << pointPosVec[iObs][2] << ")" << std::endl;
    outFile << std::setw(20) << std::setfill(' ') << "Frequency(Hz)";

    outFile << std::setw(20) << std::setfill(' ') << unit;
    outFile << std::endl;
    outFile.unsetf(std::ios::floatfield);
    outFile.precision(10);
    for (size_t iFreq = 0; iFreq < freqVec[iObs].size(); ++iFreq)
    {
      outFile << std::setw(20) << std::setfill(' ') << freqVec[iObs][iFreq];
      outFile << std::setw(20) << std::setfill(' ') << splOrPsdVec[iObs][iFreq];
      outFile << std::endl;
    }
  }
}
