library(tidyverse)
library(nlstools)
library(minpack.lm)
library(readxl)
library(broom)
library(showtext)
library(magick)
library(ggpubr)
library(lemon)


font_files() |> as_tibble() |> filter(str_detect(ps_name, "jp"))
font_add("notosans", 
         regular = "NotoSansCJKjp-Regular.otf")
theme_pubr(base_size = 12,
           base_family = "notosans") |> theme_set()
showtext_auto()

# データの読み込み
 
fname = "~/Lab_Data/学生実験_tateishi/学生実験データ/greg先生実験データ_5.12_6班_済.xlsx"
sheet = excel_sheets(fname)
range = c("A1:F85", "A1:D25", "A1:F2")
oxygen = read_xlsx(fname, 
                   sheet = sheet[1], 
                   range = range[1])

light    = read_xlsx(fname, 
                     sheet = sheet[2], 
                     range = range[2])

gww  = read_xlsx(fname, 
                 sheet = sheet[3], 
                 range = range[3])## 変数名を修正する

# tibble の処理
 
oxygen = oxygen |> 
  rename(han = "班",
         sample      = matches("サンプル"),
         net         = matches("光"),
         temperature = matches("水温"),
         min         = matches("min"),
         mgl         = matches("酸素"))

oxygen = oxygen |> drop_na()

light = light |> 
  rename(han    = "班",
         net    = "光環境",
         sample = "測定番号",
         light  = "光量子量")

light = light |> fill(han, net)

gww = gww |> 
  rename(han    = "班",
         sample = matches("サンプル"),
         gww    = matches("湿"),
         vol    = matches("酸素"))

## 光環境における光量子量の記述統計量

light = light |> 
  group_by(net) |> 
  summarise(light = mean(light))

## 光環境データと海藻データを結合する。

light = light |> 
  add_row(light = 0, net = "アルミホイル") |> 
  add_column(gww=gww$gww) |> 
  add_column(vol=gww$vol) |> 
  add_column(species = gww$海藻)

## 回帰曲線を当てはめる

pecurve = function(pmax, alpha, rd, i) {
  pmax * (1 - exp(-alpha / pmax * i)) - rd
}

oxygen = oxygen |>
  group_nest(net) |> 
  mutate(mout = map(data, ~ lm(mgl ~ min, data = .x))) |> 
  mutate(slope = map_dbl(mout, ~ coefficients(.x)[2]))

## `light` と `oxygen` を結合

alldata = full_join(light, oxygen, by = "net")

## 光合成速度を求める

alldata = alldata |> mutate(rate = slope / gww * vol)

## 光合成光曲線を当てはめる

startvalues = list(pmax = 20, alpha = 0.1, rd = 1)
preview(rate ~ pecurve(pmax, alpha, rd, light), 
        data = alldata, variable = 2,
        start = startvalues)

mfit = nlsLM(rate ~ pecurve(pmax, alpha, rd, light), data = alldata, start = startvalues)


## 非線形モデルの統計量
summary(mfit)


## Figures

cleanlabel = function(x) {sprintf("%0.1f", as.numeric(x))}

xlabel = "時間~(min)"
ylabel = "'溶存酸素濃度'~(mg~l^{-1})"

plotdata = alldata |> select(light, data) |> unnest(data) 

ggplot(plotdata) + 
  geom_point(aes(x = min, y = mgl)) + 
  geom_smooth(aes(x = min, y = mgl),
              method = "lm",
              formula = y~x) + 
  scale_x_continuous(parse(text = xlabel)) +
  scale_y_continuous(parse(text = ylabel)) +
  facet_rep_wrap(vars(light), 
             labeller = as_labeller(cleanlabel),
             ncol = 2)

pdfname = "figure-01-ogonori.pdf"
pngname = str_replace(pdfname, "pdf", "png")

w = 100
h = 150
ggsave(pdfname, width = w, height = h, units = "mm")
image_read_pdf(pdfname, density = 600) |> image_write(pngname)







newdata = alldata |> expand(light = seq(0, max(light), length = 21)) 
newdata = newdata |> mutate(fit = predict(mfit, newdata = newdata))
xlabel = "光量子量~(mu*mol~photons~m^{-2}~s^{-1})"
ylabel = "純光合成速度~(mu*g~g[{ww}]~min^{-1})"

#| fig.cap: オゴノリの光合成光曲線
ggplot() + 
  geom_point(aes(x = light, y = rate), data = alldata) +
  geom_line(aes(x = light, y = fit), data = newdata) + 
  scale_x_continuous(parse(text = xlabel)) +
  scale_y_continuous(parse(text = ylabel))


pdfname = "figure-02-ogonori.pdf"
pngname = str_replace(pdfname, "pdf", "png")

w = 100
h = 100
ggsave(pdfname, width = w, height = h, units = "mm")
image_read_pdf(pdfname, density = 600) |> image_write(pngname)






## Tables

check_p = function(x) {
  ifelse(x < 0.0001, "< 0.0001",  sprintf("%0.4f", x))
}

modelsummary = alldata |> select(light, mout) |> 
  mutate(out = map(mout, broom::glance)) |> 
  unnest(out) |> 
  select(light, adj.r.squared, statistic, df, df.residual, p.value) |> 
  unite(df, df, df.residual, sep = ", ") |> 
  arrange(light) |> 
  mutate(p.value = check_p(p.value)) 

pecurve_parameters = mfit |> coef() |> 
  as_tibble_row() |> 
  mutate(ik = pmax / alpha,
         ic = pmax / alpha * log(pmax / (pmax -rd))) |> 
  pivot_longer(everything()) 

modelsummary
pecurve_parameters

write_excel_csv(modelsummary, "ogonori-model.csv")
write_excel_csv(pecurve_parameters, "ogonori-parameters.csv")

