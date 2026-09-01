# PIEZO1 骨关节炎研究可复现代码库

本仓库 `v2.0.0` 对应最终稿的完整计算证据链。权威分析脚本、冻结结果、图表源数据、
最终稿中使用的图像、软件版本和数值校验均已分开保存。

最简执行顺序：

```text
1. 阅读 data/README.md，并取得第三方输入文件
2. python scripts/00_download_geo.py --include-large
3. python scripts/00_preflight.py --stage transcriptomic
4. python scripts/run_pipeline.py --stage transcriptomic --rscript Rscript
5. 放置 eQTLGen/FinnGen 输入后运行 genetic、validation、figures
6. python scripts/99_validate_outputs.py
```

`reference/` 是论文使用的冻结基准，重新运行生成的文件写入 `results/`、`figures/`
和 `metadata/`，不会覆盖冻结基准。遗传学第三方数据未被重新打包或重新许可；文件名、
来源、大小和 SHA-256 均记录在 `data/input_manifest.tsv`。

上传 GitHub 与 Zenodo 的中文步骤见 `docs/UPLOAD_AND_DOI_zh-CN.md`。

