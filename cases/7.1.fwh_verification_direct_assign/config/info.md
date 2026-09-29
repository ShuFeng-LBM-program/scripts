### 需要在json的"FWH"字段中增加：
* `"SourceFrequency" : "0.0017"`
* `"SourceAmplitude" : "1.0"`
格子单位。


### 需要手动修改代码编译：
`SolverLaminar::RunFwhSolvers`函数中将`fwhSolvers[grp]->update_monopole`或`fwhSolvers[grp]->update_dipole`打开，关闭`fwhSolvers[grp]->update`，即：

赋单极子流场：
```cpp
void SolverLaminar::RunFwhSolvers(const int lev)
{
    int minStep = int(flow[lev].time() / flow[lev].dt()) - 1;
    for (int grp = 0; grp < fwhSolvers.size(); ++grp)
    {
        //! for check
        const auto fwhInfoVec = Config.Fwh().get_fwhInfoVec();
        fwhSolvers[grp]->update_monopole(flow[lev].Mcr, minStep, fwhInfoVec[grp].srcFrequency, fwhInfoVec[grp].srcAmplitude);
        // fwhSolvers[grp]->update_dipole(flow[lev].Mcr, minStep, fwhInfoVec[grp].srcFrequency, fwhInfoVec[grp].srcAmplitude);
        // fwhSolvers[grp]->update(flow[lev].Mcr, minStep);
        fwhSolvers[grp]->sample(flow[lev].Mcr, minStep);
    }
}
```

赋偶极子流场：
```cpp
void SolverLaminar::RunFwhSolvers(const int lev)
{
    int minStep = int(flow[lev].time() / flow[lev].dt()) - 1;
    for (int grp = 0; grp < fwhSolvers.size(); ++grp)
    {
        //! for check
        const auto fwhInfoVec = Config.Fwh().get_fwhInfoVec();
        // fwhSolvers[grp]->update_monopole(flow[lev].Mcr, minStep, fwhInfoVec[grp].srcFrequency, fwhInfoVec[grp].srcAmplitude);
        fwhSolvers[grp]->update_dipole(flow[lev].Mcr, minStep, fwhInfoVec[grp].srcFrequency, fwhInfoVec[grp].srcAmplitude);
        // fwhSolvers[grp]->update(flow[lev].Mcr, minStep);
        fwhSolvers[grp]->sample(flow[lev].Mcr, minStep);
    }
}
```

