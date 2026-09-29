#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include <complex>
#include <iostream>
#include <functional>
#include <assert.h>

/** @brief 由时间轴计算频率轴
 * @param[in] timeVec 时间轴
 * @param[in] dataLength 2^整数(大于signalLength，dataLength/2 < signalLength)
 * @return std::vector<double> 频率轴（2^整数，小于signalLength）
 */
inline std::vector<double> cal_frequency(
  const std::vector<double>& timeVec,
  const size_t& dataLength)
{
  double startTime = timeVec[0];
  double endTime = timeVec.back();
  auto signalLength = timeVec.size();
  double dfs = double(signalLength - 1.0) / (endTime - startTime);
  std::vector<double> freqVec(dataLength / 2);
  for (size_t iFreq = 0; iFreq < dataLength / 2; ++iFreq)
  {
    freqVec[iFreq] = double(iFreq + 1.0) * dfs / double(dataLength);
    // freqVec[iFreq] = double(iFreq + 1.0) * (1.0 / (endTime - startTime)); //! 如果不考虑补零后的长度
  }
  return freqVec;
}

/**
 * @param[out] correctFactor 修正系数，弥补加窗造成的衰减：大小位2的数组指针，分别为幅值系数和能量系数
 * (source: https://community.sw.siemens.com/s/article/window-correction-factors)
 * @param[in] signalLength 信号长度
 * @param[in] windowFuncType 1: RectangularWindow; 2: HanningWindow; 3: HammingWindow; 4: KasierBesselWindow; 5: BlackmanWindow; 6: FlatTopWindow
 * @return std::vector<double> 大小为signalLength的vector，窗函数的值
 */
inline std::vector<double> cal_windowFunc_value(
  double* correctFactor,
  const size_t& signalLength,
  const int& windowFuncType)
{
  double pi = acos(-1.0);
  std::vector<double> windowFuncValue(signalLength);
  switch (windowFuncType)
  {
  case 1: // Rectangular window
    correctFactor[0] = 1.0;
    correctFactor[1] = 1.0;
    for (size_t i = 0; i < signalLength; ++i)
      windowFuncValue[i] = 1.0;
    break;
  case 2: // Hanning window
    correctFactor[0] = 2.0;
    correctFactor[1] = 1.63;
    for (size_t i = 0; i < signalLength; ++i)
      windowFuncValue[i] = 0.5 - 0.5 * cos(2.0 * pi * (i + 1.0) / signalLength);
    break;
  case 3: // Hamming window
    correctFactor[0] = 1.85;
    correctFactor[1] = 1.59;
    for (size_t i = 0; i < signalLength; ++i)
      windowFuncValue[i] = 0.54 - 0.46 * cos(2.0 * pi * (i + 1.0) / signalLength);
    break;
  case 4: // Kaiser Bessel window
    correctFactor[0] = 2.49;
    correctFactor[1] = 1.86;
    for (size_t i = 0; i < signalLength; ++i)
    {
      windowFuncValue[i] = 0.402 - 0.498 * cos(2.0 * pi * (i + 1.0) / signalLength) +
                           0.098 * cos(4.0 * pi * (i + 1.0) / signalLength) -
                           0.00123 * cos(6.0 * pi * (i + 1.0) / signalLength);
    }
    break;
  case 5: // Blackman window
    correctFactor[0] = 2.80;
    correctFactor[1] = 1.97;
    for (size_t i = 0; i < signalLength; ++i)
    {
      windowFuncValue[i] = 0.42 - 0.5 * cos(2.0 * pi * (i + 1.0) / signalLength) + 0.08 * cos(4.0 * pi * (i + 1.0) / signalLength);
    }
    break;
  case 6: // Flat Top window
    correctFactor[0] = 4.18;
    correctFactor[1] = 2.26;
    for (size_t i = 0; i < signalLength; ++i)
    {
      windowFuncValue[i] = 0.2395 - 0.448 * cos(2.0 * pi * (i + 1.0) / signalLength) +
                           0.2585 * cos(4.0 * pi * (i + 1.0) / signalLength) -
                           0.0439 * cos(6.0 * pi * (i + 1.0) / signalLength);
    }
    break;
  default:
    break;
  }
  return windowFuncValue;
}

/**
 * @param[in] windowedData 某点各时间加窗后的声压
 * windowedData.size() = dataLength (2^整数)
 */
inline void bitrp(std::vector<std::complex<double>>& windowedData)
{
  int dataLength = static_cast<int>(windowedData.size());
  int p = 0;
  for (int i = 1; i < dataLength; i *= 2)
  {
    p++;
  }
  for (int i = 0; i < dataLength; i++)
  {
    int a = i;
    int b = 0;
    for (int j = 0; j < p; j++)
    {
      b = (b << 1) + (a & 1);
      a >>= 1;
    }
    if (b > i)
    {
      std::swap(windowedData[i], windowedData[b]);
    }
  }
}

/**
 * @param[in] windowedData 某点各时间加窗后的声压
 * windowedData.size() = dataLength (2^整数次方)
 */
inline void fft(std::vector<std::complex<double>>& windowedData)
{
  double pi = acos(-1.0);
  long long dataLength = static_cast<long long>(windowedData.size());
  std::vector<std::complex<double>> w(dataLength / 2);

  bitrp(windowedData);

  double arg = -2 * pi / dataLength;
  std::complex<double> trans(cos(arg), sin(arg));

  w[0] = std::complex<double>(1.0, 0.0);
  for (int j = 1; j < dataLength / 2; j++)
  {
    w[j] = w[j - 1] * trans;
  }

  for (int m = 2; m <= dataLength; m *= 2)
  {
    for (int k = 0; k < dataLength; k += m)
    {
      for (int j = 0; j < m / 2; j++)
      {
        int index1 = k + j;
        int index2 = index1 + m / 2;
        long long t = dataLength * j / m;

        std::complex<double> trans2 = w[t] * windowedData[index2];
        std::complex<double> tmp = windowedData[index1];
        windowedData[index1] = tmp + trans2;
        windowedData[index2] = tmp - trans2;
      }
    }
  }
}

/** @brief 由时域数据计算一定范围内的频域数据
 * @param[out] freqVec 频率轴, obsNum dataLength
 * @param[out] dataFreqDomain 频域数据, obsNum dataLength
 * @param[in] timeVec 时间轴 obsNum signalLength
 * @param[in] dataTimeDomain 时域数据, obsNum signalLength
 * @param[in] windowFunc 窗函数类型 1: RectangularWindow; 2: HanningWindow; 3: HammingWindow; 4: KaiserBesselWindow; 5: BlackmanWindow; 6: FlatTopWindow
 */
inline void cal_fft(
  std::vector<std::vector<double>>& freqVec,
  std::vector<std::vector<std::complex<double>>>& dataFreqDomain,
  const std::vector<std::vector<double>>& timeVec,
  const std::vector<std::vector<double>>& dataTimeDomain,
  const int& windowFuncType)
{
  size_t obsNum = timeVec.size();
  freqVec.resize(obsNum);
  dataFreqDomain.resize(obsNum);
  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    size_t signalLength = timeVec[iObs].size();
    size_t power = 0;
    while (pow(2, power) < signalLength)
    {
      ++power;
    }
    size_t dataLength = static_cast<size_t>(pow(2, power));
    freqVec[iObs] = cal_frequency(timeVec[iObs], dataLength);

    double correctionFactor[2];
    std::vector<double> windowFuncValue = cal_windowFunc_value(correctionFactor, signalLength, windowFuncType); //! 先加窗，后补零，应该用signalLength而不是dataLength

    dataFreqDomain[iObs] = std::vector<std::complex<double>>(dataLength / 2);

    std::vector<std::complex<double>> dataFreqDomainObs(dataLength);
    for (size_t iData = 0; iData < dataLength; ++iData)
    {
      if (iData < signalLength)
        dataFreqDomainObs[iData] = std::complex<double>(windowFuncValue[iData] * dataTimeDomain[iObs][iData], 0.0);
      else
        dataFreqDomainObs[iData] = std::complex<double>(0.0, 0.0);
    }
    fft(dataFreqDomainObs);

    // 频域结果已经去掉了第一个频点为0的结果和取了前半段对称的结果
    // 单边幅值谱：|X| / (signalLength/2) 对矩形窗正确；加窗后需乘以窗函数的
    // 幅值修正系数 correctionFactor[0] = signalLength / sum(window)（即相干增益的倒数），
    // 以补偿加窗造成的幅值衰减（Hanning: ×2.0, Hamming: ×1.85, ...）。
    double scale = signalLength / 2.0;
    double ampCorrection = correctionFactor[0];
    for (size_t iFreq = 0; iFreq < dataLength / 2; ++iFreq)
    {
      dataFreqDomain[iObs][iFreq] = dataFreqDomainObs[iFreq + 1] / scale * ampCorrection;
    }
  }
}

/** @brief A计权：计算某个频率对应需要减去的声压级大小
 * @param[in] frequency 频率
 * @return double 减去的声压级大小
 */
template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
inline double cal_A_weight(T frequency)
{
  assert(frequency >= 0);
  if (abs(frequency) < 1e-15) // 可改为TF_DOUBLE_ERR
  {
    frequency = 1e-15;
  }
  double x = log10(frequency);
  double ans = (-0.190745851) * std::pow(x, 7) +
               (3.37961029) * std::pow(x, 6) +
               (-24.9512125) * std::pow(x, 5) +
               (98.6217921) * std::pow(x, 4) +
               (-221.019815) * std::pow(x, 3) +
               (257.087919) * std::pow(x, 2) +
               (-67.3872291) * std::pow(x, 1) +
               (-115.926816);

  return ans;
}

/** @brief 计算某个测点的总声压级
 * @param[in] splFreqDomainObs 该测点的频域声压级
 * @return double 该测点总声压级
 */
inline double cal_overall_soundPressureLevel(
  const std::vector<double>& splFreqDomainObs)
{
  double sumPDivPRefSqr = 0.0;
  for (size_t iFreq = 0; iFreq < splFreqDomainObs.size(); ++iFreq)
  {
    double pDivPRef = pow(10.0, splFreqDomainObs[iFreq] / 20.0);
    sumPDivPRefSqr += pDivPRef * pDivPRef;
  }
  return 20.0 * log10(sqrt(sumPDivPRefSqr));
}

/** @brief 计算三分之一倍频程对应声压级
 * @param[out] freqVecShort 1/3倍频程频率（ISO266，16Hz~16000Hz）
 * @param[out] freqVecShort 1/3倍频程对应声压级
 * @param[in] freqVec 原频率轴
 * @param[in] splVec 原声压级
 */
inline void one_third_octave_band(
  std::vector<double>& freqVecShort,
  std::vector<double>& splVecShort,
  const std::vector<double>& freqVec,
  const std::vector<double>& splVec)
{
  std::vector<double> octaves = {16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0,
                                 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0,
                                 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0,
                                 4000.0, 5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0};
  size_t dataLength = freqVec.size();
  size_t dataLengthShort = octaves.size();
  freqVecShort = octaves;
  splVecShort.assign(dataLengthShort, 0.0);
  for (size_t oct = 0; oct < dataLengthShort; ++oct)
  {
    double freqUpper = octaves[oct] * sqrt(pow(2.0, 1.0 / 3.0));
    double freqLower = octaves[oct] / sqrt(pow(2.0, 1.0 / 3.0));
    double sumPower = 0.0;
    for (size_t iFreq = 0; iFreq < dataLength; ++iFreq)
    {
      if (freqVec[iFreq] < freqLower)
        continue;
      else if (freqVec[iFreq] >= freqLower && freqVec[iFreq] <= freqUpper)
      {
        sumPower += pow(10.0, splVec[iFreq] / 10.0);
      }
      else if (freqVec[iFreq] > freqUpper)
        break;
    }
    splVecShort[oct] = 10.0 * log10(sumPower);
  }
}

/** @brief 取出频率范围内的频率轴和对应频域数据
 * @param[out] freqShort 指定范围频率轴
 * @param[out] dataFreqDomainShort 指定频率范围的频域数据
 * @param[in] freq 原频率轴
 * @param[in] dataFreqDomain 原频域数据
 * @param[in] minFreq 最小频率
 * @param[in] maxFreq 最大频率
 */
inline void apply_freq_limits(
  std::vector<double>& freqShort,
  std::vector<double>& dataFreqDomainShort,
  const std::vector<double>& freq,
  const std::vector<double>& dataFreqDomain,
  const double& minFreq,
  const double& maxFreq)
{
  // apply minFreq, maxFreq limits
  size_t dataLength = freq.size();
  size_t minFreqIdx = 0;
  size_t maxFreqIdx = 0;
  size_t iFreq = 0;
  while (freq[iFreq] < minFreq && iFreq < dataLength)
  {
    minFreqIdx++;
    iFreq++;
  }
  iFreq = 0;
  while (freq[iFreq] < maxFreq && iFreq < dataLength)
  {
    maxFreqIdx++;
    iFreq++;
  }
  std::vector<double>::const_iterator firstFreqIt = freq.begin() + minFreqIdx;
  std::vector<double>::const_iterator lastFreqIt = freq.begin() + maxFreqIdx + 1;
  freqShort.assign(firstFreqIt, lastFreqIt);

  std::vector<double>::const_iterator firstDataIt = dataFreqDomain.begin() + minFreqIdx;
  std::vector<double>::const_iterator lastDataIt = dataFreqDomain.begin() + maxFreqIdx + 1;
  dataFreqDomainShort.assign(firstDataIt, lastDataIt);
}

/** @brief 计算两个三维向量的点乘
 * @param[in] p1 向量p1
 * @param[in] p2 向量p2
 * @return double 向量p1和向量p2点乘的结果
 */
inline double dotproduct(
  const double* p1,
  const double* p2,
  const size_t& length)
{
  double ans = 0.0;
  for (size_t i = 0; i < length; ++i)
  {
    ans += p1[i] * p2[i];
  }
  return ans;
}

/** @brief 计算两个数组的卷积
 * @param[in] arrayA 数组A
 * @param[in] arrayB 数组B
 * @param[in] mode "full", "same", "valid" (only "valid" mode is supported here)
 * @return std::vector<double> 数组A和数组B的卷积
 */
inline std::vector<double> convolve(
  const std::vector<double>& arrayA,
  const std::vector<double>& arrayB,
  const std::string& mode = "valid")
{
  std::vector<double> arrShorter = arrayA;
  std::vector<double> arrLonger = arrayB;
  if (arrayA.size() > arrayB.size())
  {
    arrShorter = arrayB;
    arrLonger = arrayA;
  }
  auto lengthShorter = arrShorter.size();
  auto lengthLonger = arrLonger.size();
  size_t length;
  std::vector<double> ans;

  if (mode == "valid")
  {
    std::reverse(arrShorter.begin(), arrShorter.end());
    length = lengthLonger - lengthShorter + 1;
    ans.resize(length);
    std::vector<double>::const_iterator itLonger = arrLonger.begin();
    std::vector<double> arrLongerShort;
    for (size_t i = 0; i < length; ++i)
    {
      arrLongerShort.assign(itLonger + i, itLonger + i + lengthShorter);
      ans[i] = dotproduct(arrShorter.data(), arrLongerShort.data(), lengthShorter);
    }
  }
  return ans;
}
/** @brief 由频域声压计算声压级，考虑倍频和计权
 * @param[out] freqWeightedBandVec 新频率轴（考虑倍频、计权）
 * @param[out] splWeightedBandVec 频域声压级（考虑倍频、计权）
 * @param[out] freqSmoothedVec 频率轴（窄带，考虑计权、光顺）
 * @param[out] splWeightedSmoothedVec 频域声压级（窄带，考虑计权、光顺）
 * @param[in] freqVec 频率轴
 * @param[in] preFreqDomainVec 频域压力
 * @param[in] minFreq 最小频率（用于计算总声压级）
 * @param[in] maxFreq 最大频率（用于计算总声压级）
 * @param[in] band 倍频程
 * @param[in] weighting 计权方式
 * @param[in] nSmooth 光顺系数（大于1，默认为5）
 * @return std::vector<double> 各测点的总声压级（窄带，考虑计权）obsNum
 */
inline std::vector<double> pressure_to_spl_smoothed(
  std::vector<std::vector<double>>& freqWeightedBandVec,
  std::vector<std::vector<double>>& splWeightedBandVec,
  std::vector<std::vector<double>>& freqSmoothedVec,
  std::vector<std::vector<double>>& splWeightedSmoothedVec,
  const std::vector<std::vector<double>>& freqVec,
  const std::vector<std::vector<std::complex<double>>>& preFreqDomainVec,
  const double& minFreq,
  const double& maxFreq,
  const std::string& band,
  const std::string& weighting,
  const int& nSmooth = 5)
{
  double refPressure = 2.0E-5; // Pa
  size_t obsNum = freqVec.size();
  freqWeightedBandVec.resize(obsNum);
  splWeightedBandVec.resize(obsNum);
  freqSmoothedVec.resize(obsNum);
  splWeightedSmoothedVec.resize(obsNum);

  std::function<double(double)> weightingFunc;
  if (weighting == "A-weight")
  {
    weightingFunc = [](double frequency)
    { return cal_A_weight(frequency); };
  }
  else if (weighting == "none")
  {
    weightingFunc = [](double frequency)
    { return 0.0; };
  }

  // calculate non-weighted/weighted SPL & overall SPL
  std::vector<double> oaSplWeightedSmoothedVec(obsNum, 0.0);
  std::vector<std::vector<double>> splWeightedSmoothedVecTmp(obsNum); // 用于计算总声压级
  std::vector<std::vector<double>> freqSmoothedVecTmp(obsNum);        // 用于计算总声压级
  std::vector<double> weights(nSmooth, 1.0 / double(nSmooth));
  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    size_t dataLength = freqVec[iObs].size();
    splWeightedBandVec[iObs].resize(dataLength);
    for (size_t iFreq = 0; iFreq < dataLength; ++iFreq)
    {
      splWeightedBandVec[iObs][iFreq] = 20.0 * log10(std::abs(preFreqDomainVec[iObs][iFreq]) / sqrt(2.0) / refPressure); //! 除以根号2，用有效幅值
      splWeightedBandVec[iObs][iFreq] += weightingFunc(freqVec[iObs][iFreq]);
    }
    freqSmoothedVec[iObs] = convolve(freqVec[iObs], weights);
    splWeightedSmoothedVec[iObs] = convolve(splWeightedBandVec[iObs], weights);
    // 考虑频率范围计算总声压级
    apply_freq_limits(freqSmoothedVecTmp[iObs], splWeightedSmoothedVecTmp[iObs], freqSmoothedVec[iObs], splWeightedSmoothedVec[iObs], minFreq, maxFreq);
    oaSplWeightedSmoothedVec[iObs] = cal_overall_soundPressureLevel(splWeightedSmoothedVecTmp[iObs]);
  }

  // consider octave band
  if (band == "narrowBand")
  {
    for (size_t iObs = 0; iObs < obsNum; ++iObs)
    {
      freqWeightedBandVec[iObs] = freqSmoothedVec[iObs];
      splWeightedBandVec[iObs] = splWeightedSmoothedVec[iObs];
    }
  }
  else if (band == "oneThirdOctaveBand")
  {
    for (size_t iObs = 0; iObs < obsNum; ++iObs)
    {
      one_third_octave_band(freqWeightedBandVec[iObs], splWeightedBandVec[iObs], freqSmoothedVec[iObs], splWeightedSmoothedVec[iObs]);
    }
  }
  return oaSplWeightedSmoothedVec;
}

/** @brief 频域数据计算功率谱密度
 * @param[in] dataFreqDomain 频域数据（已经除以(signalLength/2)）
 * @param[in] signalLength 时域数据长度
 * @param[in] dt 时域时间步长
 * @return std::vector<double> 功率谱密度
 */
inline std::vector<double> psd_raw(
  const std::vector<std::complex<double>>& dataFreqDomain,
  const size_t& signalLength,
  const double& dt)
{
  std::vector<double> psdRaw;
  for (size_t iData = 0; iData < dataFreqDomain.size(); ++iData)
  {
    double abs = std::abs(dataFreqDomain[iData]) / sqrt(2);
    psdRaw.push_back(abs * abs * signalLength * dt / 2.0);
  }
  return psdRaw;
}

/** @brief 由频域数据计算光顺后的功率谱密度
 * @param[out] freqSmoothed 光顺后的频率轴
 * @param[out] psdSmoothed 光顺后的功率谱密度
 * @param[in] freq 频率轴
 * @param[in] dataFreqDomain 频域数据（已经除以(signalLength/2)）
 * @param[in] signalLength 时域数据长度
 * @param[in] dt 时域时间步长
 * @param[in] nSmooth 光顺系数（大于1，默认为5）
 */
inline void psd_smoothed(
  std::vector<double>& freqSmoothed,
  std::vector<double>& psdSmoothed,
  const std::vector<double>& freq,
  const std::vector<std::complex<double>>& dataFreqDomain,
  const size_t& signalLength,
  const double& dt,
  const int& nSmooth = 5)
{
  auto psdRaw = psd_raw(dataFreqDomain, signalLength, dt);
  std::vector<double> weights(nSmooth, 1.0 / double(nSmooth));
  psdSmoothed = convolve(psdRaw, weights);
  freqSmoothed = convolve(freq, weights);
}

/** @brief 由频域声压计算功率谱密度级并光顺处理
 * @param[out] freqSmoothedVec  新频率轴
 * @param[out] psdSmoothedVec 光顺后的功率谱密度级
 * @param[in] freqVec 原频率轴
 *
 */
inline void pressure_to_psd_smoothed(
  std::vector<std::vector<double>>& freqSmoothedVec,
  std::vector<std::vector<double>>& psdSmoothedVec,
  const std::vector<std::vector<double>>& freqVec,
  const std::vector<std::vector<std::complex<double>>>& dataFreqDomainVec,
  const std::vector<std::vector<double>>& timeVec,
  const double& minFreq,
  const double& maxFreq,
  const int& nSmooth = 5)
{
  const size_t obsNum = freqVec.size();
  freqSmoothedVec.resize(obsNum);
  psdSmoothedVec.resize(obsNum);

  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    size_t signalLength = timeVec[iObs].size();
    double dt = (timeVec[iObs].back() - timeVec[iObs][0]) / double(signalLength - 1);
    psd_smoothed(freqSmoothedVec[iObs], psdSmoothedVec[iObs], freqVec[iObs], dataFreqDomainVec[iObs], signalLength, dt, nSmooth);
  }
}

/** @brief 由计权的声压级计算功率谱密度(级)并光顺处理
 * @param[out] psdWeightedSmoothedVec 光顺后的功率谱密度
 * @param[out] psdWeightedSmoothedLevelVec 光顺后的功率谱密度级
 * @param[in] freqSmoothedVec 频率轴
 * @param[in] splWeightedVec 计权的声压级
 */
inline void spl_to_psd_smoothed(
  std::vector<std::vector<double>>& psdWeightedSmoothedVec,
  std::vector<std::vector<double>>& psdWeightedSmoothedLevelVec,
  const std::vector<std::vector<double>>& freqSmoothedVec,
  const std::vector<std::vector<double>>& splWeightedSmoothedVec)
{
  size_t obsNum = freqSmoothedVec.size();
  psdWeightedSmoothedVec.resize(obsNum);
  psdWeightedSmoothedLevelVec.resize(obsNum);

  for (size_t iObs = 0; iObs < obsNum; ++iObs)
  {
    auto dataLengthSM = freqSmoothedVec[iObs].size();
    double dfSM = (freqSmoothedVec[iObs].back() - freqSmoothedVec[iObs][0]) / double(dataLengthSM - 1);
    for (size_t iFreq = 0; iFreq < dataLengthSM; ++iFreq)
    {
      double psdWeightedSmoothedLvl = splWeightedSmoothedVec[iObs][iFreq] - 10.0 * log10(dfSM);
      psdWeightedSmoothedLevelVec[iObs].push_back(psdWeightedSmoothedLvl);
      psdWeightedSmoothedVec[iObs].push_back(pow(10.0, psdWeightedSmoothedLvl / 10.0));
    }
  }
}
