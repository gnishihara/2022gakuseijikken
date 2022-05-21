# オゴノリの光合成光曲線の解析
# 2022-05-21


# パッケージ読み込み

library(tidyverse)
library(nlstools)
library(ggpubr)
library(magick)
library(lemon)
library(readxl)
library(minpack.lm)

library(broom)
library(showtext)

# パソコンにインストールされたフォントをさがす
font.files() |> as_tibble() |> 
  filter(str_detect(ps_name, "jp"))

font_add(family = "forplot",
         regular = "NotoSansCJKjp-Regular.otf")
showtext_auto()

# データの読み込み

fname = "Data/greg先生実験データ_5.12_6班_済.xlsx"

sheets = excel_sheets(fname)
range = c("A1:F85",
          "A1:D26",
          "A1:F2")

# いちどに全データを読み込む場合
# df = tibble(fname, sheets, range)
# df |> mutate(data = pmap(list(fname, sheets, range), read_xlsx))

#　ここのシートを読み込む場合
oxygen  = read_xlsx(path = fname,  sheet = sheets[1], range = range[1])
light   = read_xlsx(path = fname,  sheet = sheets[2], range = range[2])
seaweed = read_xlsx(path = fname,  sheet = sheets[3], range = range[3])

# 列名を英語に変換する


oxygen = oxygen |> 
  rename(han    = "班",
         sample = "サンプル",
         min    = matches("時間"),
         mgl    = `酸素（mg/L)`,
         temperature = matches("水温"),
         light = "光環境")

light = light |> 
  rename(han = "班",
         light = "光環境",
         ppfd = "光量子量")

seaweed = seaweed |> 
  rename(han = "班",
         kaisou = "海藻",
         sample = "サンプル",
         gww = matches("湿重量"),
         vol = matches("容量"))


# NA を fill する

light = light |> fill(han, light)

light |> print(n = Inf)

# light の平均PPFDを求める。

light = light |> 
  group_by(light) |> 
  summarise(ppfd_mean = mean(ppfd))


oxygen
light
seaweed

tmp = bind_cols(light, seaweed)


alldata = left_join(oxygen, tmp, by = c("light", "han", "sample"))

alldata |> print(n = Inf)

# alldata から NA を外す

alldata = alldata |> drop_na()

alldata |> print(n = Inf)


# 光合成速度の計算

alldata =   alldata |> group_nest(light, ppfd_mean, gww, vol)

# 自作関数

fit_lm = function(data) {
  lm(mgl ~ min, data = data)
}

fit_lm

get_slope = function(model) {
  coefficients(model)[2]
}

get_slope

# 自作関数をつかいます。

alldata = alldata |> mutate(model = map(data, fit_lm))

alldata = alldata |> mutate(slope = map(model, get_slope))

alldata = alldata |> unnest(slope)

# map() でなにをしているのか？ ###################################

light0 = alldata |> slice(1) |> unnest(data)
model1 = fit_lm(light0)
slope1 = get_slope(model)

light2 = alldata |> slice(2) |> unnest(data)
model2 = fit_lm(light2)
slope2 = get_slope(model2)
##################################################################

# 光合成光速度の統計量
# adj.r.squared: 調整済み決定係数
# statistic: F値
# df, df.residual: F値の自由度
alldata |> 
  mutate(stats = map(model, glance)) |> 
  unnest(stats) |> 
  select(light, ppfd_mean, adj.r.squared, statistic, df, df.residual) |> 
  mutate(pvalue = df(statistic, df, df.residual))

# 光合成速度の回帰曲線の図



alldata_fig = alldata |> unnest(data)

xlabel = "Time~(min)"
ylabel = "Dissolved~oxygen~(mg~l^{-1})"
clabel = "Mean~PPFD"

fix_legend = function(x) {
  as.numeric(x) |> format(digits = 1)
}

ggplot(alldata_fig) + 
  geom_point(aes(x = min, y = mgl, 
                 color = as.factor(ppfd_mean))) +
  geom_smooth(aes(x = min, y = mgl,
                  color = as.factor(ppfd_mean)),
              method = "lm",
              formula = y ~ x, 
              se = FALSE) +
  scale_x_continuous(parse(text = xlabel)) +
  scale_y_continuous(parse(text = ylabel)) +
  scale_color_viridis_d(parse(text = clabel),
                        end = 0.9,
                        labels = fix_legend) +
  labs(title = "オゴノリの酸素濃度") +
  theme_pubr(base_family = "forplot")

pdfname = "ogonori_fig01.pdf"  
pngname = "ogonori_fig01.png"  

ggsave(pdfname, width = 150, height = 100, units ="mm")
image_read_pdf(pdfname) |> image_write(pngname)


## 光合成光曲線


pecurve = function(pmax, alpha, rd, ppfd) {
  pmax * (1-exp(-alpha / pmax * ppfd)) - rd
}

alldata = alldata |> 
  mutate(rate = slope * vol / gww)


start = list(pmax = 20, rd = 2, alpha = 0.1)

preview(rate ~ pecurve(pmax, alpha, rd, ppfd_mean),
        data = alldata,
        start = start,
        variable = 2)

model_ogo = nlsLM(rate ~ pecurve(pmax, alpha, rd, ppfd_mean),
                  data = alldata,
                  start = start)

summary(model_ogo)

# 光飽和点
ik_fn = function(alpha, pmax) {
  x = pmax / alpha 
  names(x) = "ik"
  x
}

# 光補償点
ic_fn = function(alpha, pmax, rd) {
  x = pmax / alpha * log(pmax / (pmax - rd))
  names(x) = "ic"
  x
}

cfs_ogo = coef(model_ogo)

ik = ik_fn(alpha = cfs_ogo["alpha"], 
           pmax = cfs_ogo["pmax"])

ic = ic_fn(alpha = cfs_ogo["alpha"], 
           pmax = cfs_ogo["pmax"],
           rd = cfs_ogo["rd"])
ic 
ik

c(cfs_ogo, ic, ik)



## オゴノリの光合成光曲線

newdata = tibble(ppfd_mean = seq(0, 500, length = 21))

newdata = newdata |> 
  mutate(fit = predict(model_ogo, newdata = newdata))

xlabel = "PPFD~(mu*mol~photons~m^{-2}~s^{-1})"
ylabel = "NP~(mu*g~g[{ww}]^{-1}~min^{-1})"

ggplot() +
  geom_point(aes(x = ppfd_mean, y = rate),
             data = alldata) +
  geom_line(aes(x = ppfd_mean, y = fit),
            data = newdata) +
  scale_x_continuous(parse(text = xlabel)) +
  scale_y_continuous(parse(text = ylabel)) +
  labs(title = "オゴノリの光合成光曲線") +
  theme_pubr(base_family = "forplot")


pdfname = "ogonori_fig02.pdf"  
pngname = "ogonori_fig02.png"  

ggsave(pdfname, width = 150, height = 100, units ="mm")
image_read_pdf(pdfname) |> image_write(pngname)








