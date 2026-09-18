# 水利工程 AI 技能库

面向水利水电工程设计的 AI 技能（Skills）集合，共 **35 个技能**，全部可在支持
Agent Skills 规范的环境（Claude Code / WorkBuddy / OpenClaw 等）中直接使用。

**收录原则**：只放能算、能用、算得对的东西。算法讲得清来源，实现讲得清依据，
错误标注在明面上——**不把错误的知识传给下一个人**。

---

## 目录结构

```
skills/
├── slcalc/          # 《水利水电工程设计计算程序集》Python 改造系列（32 个）
└── original/        # 独立原创技能（3 个）
```

每个技能目录内均含 `SKILL.md`（元数据 + 用法 + 算法说明）与配套实现。

---

## 一、《水利水电工程设计计算程序集》改造系列（32 个）

原程序集由**乌鲁木齐正海水利科技有限公司**开发（张校正教授级高工技术总负责），
2025 年 3 月公之于众。该程序集是国家《水工设计手册》选定的"水工设计常用软件"之一，
1990 年获能源部水利部水电规划设计总院"七五"计算机应用优秀成果一等奖，
2010 年获国家版权局软件著作权登记。

本系列为**现代环境下的 Python 改造实现**，覆盖工程水文计算、水利水能计算、
工程水力学计算等专业方向。

| 技能 | slug | 版本 |
|---|---|---|
| [重力式挡土墙计算程序（G-6）](skills/slcalc/slcalc-dangtuqiang/) | `slcalc-dangtuqiang` | v1.0.1 |
| [分析经验单位线及汇流计算程序（A-7）](skills/slcalc/slcalc-fenxi-jingyan-danweixian-ji-huiliu-jisuan/) | `slcalc-fenxi-jingyan-danweixian-ji-huiliu-jisuan` | v1.0.0 |
| [分析马司京根法演算参数程序（A-8）](skills/slcalc/slcalc-fenxi-masijinggenfa-yansuan-canshu/) | `slcalc-fenxi-masijinggenfa-yansuan-canshu` | v1.0.0 |
| [美国气象局溃坝洪水预报简化模型程序（C-1）](skills/slcalc/slcalc-kuiba-hongshui-yubao-jianhua-moxing/) | `slcalc-kuiba-hongshui-yubao-jianhua-moxing` | v1.0.0 |
| [两种供水保证的水库径流调节计算程序（C-3）](skills/slcalc/slcalc-liangzhong-gongshui-baozheng-shuiku-jingliu-tiaojie-jisuan/) | `slcalc-liangzhong-gongshui-baozheng-shuiku-jingliu-tiaojie-jisuan` | v1.0.0 |
| [马斯京根模型最优参数估计程序（A-14）](skills/slcalc/slcalc-masijinggen-moxing-zuiyou-canshu-guiji/) | `slcalc-masijinggen-moxing-zuiyou-canshu-guiji` | v1.0.0 |
| [马司京根法分段连续演算程序（A-9）](skills/slcalc/slcalc-masijinggenfa-fenduan-lianxu-yansuan/) | `slcalc-masijinggenfa-fenduan-lianxu-yansuan` | v1.0.0 |
| [常用断面渠道水力学计算程序（D-6）](skills/slcalc/slcalc-qudao-shuilixue/) | `slcalc-qudao-shuilixue` | v1.0.0 |
| [三参数幂函数曲线拟合程序（A-2）](skills/slcalc/slcalc-sancanshu-mihanshu-quxiannihe/) | `slcalc-sancanshu-mihanshu-quxiannihe` | v1.0.0 |
| [水电站等出力调节计算(试算法)程序（C-7）](skills/slcalc/slcalc-shuidianzhan-dengchuli-tiaojie-shisuanfa/) | `slcalc-shuidianzhan-dengchuli-tiaojie-shisuanfa` | v1.0.0 |
| [水电站等出力调节计算(图解法)程序（C-9）](skills/slcalc/slcalc-shuidianzhan-dengchuli-tiaojie-tujiefa/) | `slcalc-shuidianzhan-dengchuli-tiaojie-tujiefa` | v1.0.0 |
| [水电站等流量调节计算通用程序（C-4）](skills/slcalc/slcalc-shuidianzhan-dengliuliang-tiaojie-jisuan/) | `slcalc-shuidianzhan-dengliuliang-tiaojie-jisuan` | v1.0.0 |
| [水电站定出力调节计算(试算法)程序（C-6）](skills/slcalc/slcalc-shuidianzhan-dingchuli-tiaojie-shisuanfa/) | `slcalc-shuidianzhan-dingchuli-tiaojie-shisuanfa` | v1.0.0 |
| [水电站定出力调节计算(图解法)程序（C-8）](skills/slcalc/slcalc-shuidianzhan-dingchuli-tiaojie-tujiefa/) | `slcalc-shuidianzhan-dingchuli-tiaojie-tujiefa` | v1.0.0 |
| [水电站和水利枢纽年调节计算程序（C-14）](skills/slcalc/slcalc-shuidianzhan-shuili-shuniu-niantiaojie-jisuan/) | `slcalc-shuidianzhan-shuili-shuniu-niantiaojie-jisuan` | v1.0.0 |
| [水电站蓄水期调节计算程序（C-5）](skills/slcalc/slcalc-shuidianzhan-xushuiqi-tiaojie-jisuan/) | `slcalc-shuidianzhan-xushuiqi-tiaojie-jisuan` | v1.0.0 |
| [水库调洪演算的数值解程序（C-2）](skills/slcalc/slcalc-shuiku-tiaohong-yansuan/) | `slcalc-shuiku-tiaohong-yansuan` | v1.0.0 |
| [水文频率计算程序（A-3）](skills/slcalc/slcalc-shuiwen-pinlv/) | `slcalc-shuiwen-pinlv` | v1.0.0 |
| [水文频率计算程序·连续/不连续系列（A-3X）](skills/slcalc/slcalc-shuiwen-pinlv-x/) | `slcalc-shuiwen-pinlv-x` | v1.0.0 |
| [水文系列代表性分析程序（A-4）](skills/slcalc/slcalc-shuiwen-xilie-daibiaoxing-fenxi/) | `slcalc-shuiwen-xilie-daibiaoxing-fenxi` | v1.0.0 |
| [水闸水力学计算程序（D-2）](skills/slcalc/slcalc-shuizha-shuilixue-jisuan/) | `slcalc-shuizha-shuilixue-jisuan` | v1.0.0 |
| [随机水文AR(P)模型分析程序（A-13）](skills/slcalc/slcalc-suiji-shuiwen-ar-moxing-fenxi/) | `slcalc-suiji-shuiwen-ar-moxing-fenxi` | v1.0.0 |
| [推求法计算天然河道水面曲线程序（D-14A）](skills/slcalc/slcalc-tianran-hedao-shuimianquxian/) | `slcalc-tianran-hedao-shuimianquxian` | v1.0.0 |
| [调水工程水利经济计算程序（C-12）](skills/slcalc/slcalc-tiaoshui-gongcheng-shuili-jingji-jisuan/) | `slcalc-tiaoshui-gongcheng-shuili-jingji-jisuan` | v1.0.0 |
| [同频率缩放设计洪水过程线程序（A-5X）](skills/slcalc/slcalc-tongpinlv-suofang-sheji-hongshui-guochengxian/) | `slcalc-tongpinlv-suofang-sheji-hongshui-guochengxian` | v1.0.0 |
| [推理公式法计算洪峰流量程序（A-10）](skills/slcalc/slcalc-tuilu-gongshifa-jisuan-hongfengliuliang/) | `slcalc-tuilu-gongshifa-jisuan-hongfengliuliang` | v1.0.0 |
| [推求流域时段平均面雨量程序（A-6）](skills/slcalc/slcalc-tuiqiu-liuyu-shiduan-pingjun-mianyuliang/) | `slcalc-tuiqiu-liuyu-shiduan-pingjun-mianyuliang` | v1.0.0 |
| [无调节水电站水能计算程序（C-11）](skills/slcalc/slcalc-wutiaojie-shuidianzhan-shuineng-jisuan/) | `slcalc-wutiaojie-shuidianzhan-shuineng-jisuan` | v1.0.0 |
| [下渗曲线产流计算程序（A-12）](skills/slcalc/slcalc-xiashen-quxian-chanliu-jisuan/) | `slcalc-xiashen-quxian-chanliu-jisuan` | v1.0.0 |
| [由直方图所限定的曲线拟合程序（A-15）](skills/slcalc/slcalc-zhifangtu-xianxian-quxian-nihe/) | `slcalc-zhifangtu-xianxian-quxian-nihe` | v1.0.0 |
| [直线相关计算程序（A-1）](skills/slcalc/slcalc-zhixian-xiangguan/) | `slcalc-zhixian-xiangguan` | v1.0.0 |
| [最大24小时洪量计算程序（A-11）](skills/slcalc/slcalc-zuida-24xiaoshi-hongliang-jisuan/) | `slcalc-zuida-24xiaoshi-hongliang-jisuan` | v1.0.0 |

---

## 二、独立原创技能（3 个）

| 技能 | slug | 版本 |
|---|---|---|
| [高原·水利工程专家](skills/original/gaoyuan-shuili-gongcheng-zhuanjia/) | `gaoyuan-shuili-gongcheng-zhuanjia` | v1.0.0 |
| [小说九维评价标准](skills/original/novel-eval-9dim/) | `novel-eval-9dim` | v1.0.0 |
| [水库双辅助曲线调洪演算](skills/original/reservoir-flood-routing/) | `reservoir-flood-routing` | v1.0.0 |

---

## 使用方法

**方式一：直接放入技能目录**

```bash
git clone https://github.com/<your-account>/shuili-skills.git
cp -r shuili-skills/skills/* <你的技能目录>/
```

各环境默认技能目录：
- Claude Code：`~/.claude/skills/`
- WorkBuddy / OpenClaw：`~/.workbuddy/skills/`

**方式二：按需取用**

每个技能目录可独立复制，`SKILL.md` 顶部 frontmatter 说明触发方式与依赖。

---

## 授权与来源

- **本仓库含两类来源不同的内容，授权条款不同，请务必阅读 [DECLARATION.md](DECLARATION.md)。**
- 改造系列（`skills/slcalc/`）：算法与原始数据版权归**原开发单位**所有，
  本仓库提供的是 Python 改造实现。
- 原创技能（`skills/original/`）：版权归作者所有。

---

## 问题反馈

发现算错、跑不通、或者哪条公式与规范对不上——请直接开
[Issue](https://github.com/<your-account>/shuili-skills/issues) 说明。

**算得对是这里唯一的标准。** 任何计算错误、与规范/手册不一致之处，
欢迎指出，会核对后修正并在提交记录中留痕。

---

## 免责声明

本仓库技能供工程技术人员学习与辅助计算使用。**所有计算结果在用于实际工程前，
必须由具备相应资质的工程技术人员按现行有效的国家标准、行业规范复核确认。**
作者不对因使用本仓库内容造成的任何后果承担责任。

---

*本仓库为技术存档与共享用途，随技能修订持续更新。*
