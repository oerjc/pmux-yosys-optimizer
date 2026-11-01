# 从源码重新构建

本入口从包内原始归档重新编译一个匹配的 Windows 双版本。它先构建公共 ABC 和原版 Yosys，再应用候选归档中经过范围审查的 PMUX 文件、两处 CMake 注册与一次 `synth_intel` 调用，构建优化版。其他 pass、ABC9 与运行时均使用上游源码。

本地基线为公开 v0.69 提交 `9f75ca1f9834a39a863915b5dae0c7b1e33533bc`；ABC 提交为 `a51bf4a7bb38ce06ea3a544b1463b40832213917`。官方说明中的 `0.69+24 d0e71cfb7-dirty` 含未公开修改，因此本包不声称与该二进制完全一致。重新编译得到的新二进制也不自动继承包内正式测试成绩。

## 外部条件

- Windows x64、Python 3.12 或更新版本。
- 已安装的 MSYS2 UCRT64：GCC/G++、Python、GNU make、Flex、Bison，以及对应 GCC 运行库。MSYS2 根路径须为不含空格的 ASCII 路径。
- 一个尚不存在的构建目录，以及一个空闲的 P 至 Z 盘符，默认 P。脚本仅撤销自己创建的 `subst` 映射。
- 充足的磁盘空间和时间。双构建远比运行现有二进制耗时；正式性能采样期间不要并行编译。

原开发机版本为 GCC/G++ 15.2.0、Bison 3.8.2、Flex 2.6.4、MSYS2 Python 3.12.11。CMake 4.4.4、Ninja 1.13.2 及所需源码归档已在 `vendor` 中；入口不会联网下载或安装依赖。依赖包与候选源码摘要在 `build-inputs.json`，许可在 `licenses`。

## 操作

在提交包根目录打开 PowerShell，按本机真实路径替换示例：

```powershell
& .\source\build_windows.ps1 -Python 'C:\Python312\python.exe' `
  -BuildRoot 'D:\pmux-rebuild-01' -MsysRoot 'D:\MSYS2' -Jobs 2 -CheckOnly
```

`-CheckOnly` 只核对归档摘要、允许的源码差异、外部工具和空闲路径，不编译、不生成成绩。通过后去掉 `-CheckOnly`，其余参数不变，开始实际构建。若 P 盘已使用，可显式增加 `-Drive Q`。失败时保留新构建目录与日志以便检查；再次构建应另选新目录。

构建期间两臂共用相同源码路径、编译器、CMake 开关、ABC 和依赖；先保存原版的编译记录，再在同一路径叠加四项已审查的源码差异。脚本检查所有公共编译命令一致，候选仅增加 `pmux_opt.cc` 编译，原生 `zinit` 在两臂均启用。配置声明了 Release 与 `-O2 -DNDEBUG`，但 v0.69 自身目标规则实际使用的优化选项以保存的 `compile_commands.json` 为准；原构建的公共 Yosys 命令为 `-O3`，不能仅凭配置参数宣称最终使用了 `-O2`。

成功输出位于：

```text
<BuildRoot>/work/final-toolchain/install-baseline/
<BuildRoot>/work/final-toolchain/install-candidate/
<BuildRoot>/build-records/baseline/
<BuildRoot>/build-records/candidate/
<BuildRoot>/rebuild-summary.json
```

运行树包含 `bin/yosys.exe`、公共 `yosys-abc.exe`、必要 DLL、`share/yosys` 数据和许可证。`rebuild-summary.json` 记录本次实际摘要、外部工具摘要与比较结果；其 `BUILT` 只表示编译和轻量运行检查成功，未进行综合、等价证明或性能测量。

## 验证新构建

仍从提交包根目录执行：

```powershell
& .\scripts\run-validation.ps1 -Python 'C:\Python312\python.exe' `
  -SuiteRoot 'D:\Tools\oss-cad-suite' -Bash 'D:\MSYS2\usr\bin\bash.exe' `
  -RebuildRoot 'D:\pmux-rebuild-01' -RunName 'rebuilt-formal-01'
```

`-RebuildRoot` 显式切换到新双版本，先核对其摘要与本次构建清单，再生成全新的九例正式计划和证明结果。若不提供此参数，验证入口使用提交包原附的两个运行树。结果写入提交包 `runs/`，不会替换 `reports/final/` 中的历史正式证据。只有实际重新通过的结果才适用于新二进制。

构建入口已经按本项目真实配置整理，但是否已在新目录完整演练必须以包内 `build/PORTABLE_BUILD_CHECK.md` 的实际记录为准；语法或准备检查不等于完成重编译。
