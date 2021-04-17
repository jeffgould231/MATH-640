btc_first <- read_delim("BTC_2021-02-04.txt", col_names = F, delim = ",")
colnames(btc_first) <- c("TimeStamp", "Price", "Volume", "Exchange")

btc_first %>%
  mutate(Time = lubridate::as_datetime(TimeStamp))


vol_candle_size <- 500

vol_candle <- btc_first %>%
       mutate(CumVol = cumsum(Volume)) %>%
  mutate(VolBar = ceiling(lag(CumVol) / vol_candle_size)) %>%
    group_by(VolBar) %>%
    summarise(Start = first(TimeStamp),
              End = last(TimeStamp),
              Open = first(Price),
              Close = last(Price),
              High = max(Price),
              Low = min(Price),
              Variance = sd(Price)^2,
              Return = log(last(Price))- log(first(Price)),
              Trades = n(),
              TotalVol = sum(Volume),
              AvgVol = mean(Volume))

min(vol_candle$TotalVol)
max(vol_candle$TotalVol)
mean(vol_candle$TotalVol)
quantile(vol_candle$TotalVol)

ggplot(vol_candle, aes(x = VolBar)) +
  geom_candlestick(aes(open = Open, close = Close, high = High, low = Low))
  geom_density(aes(x = Trades))


View(tq_get("AAPL", from = "2013-01-01", to = "2016-12-31"))

library(jsonlite)
library(glue)

# create a function to retrieve daily data
retreive_daily_data <- function(pair, filename) {
  url = glue("https://api.pro.coinbase.com/products/{pair}/candles?granularity=86400")
  columnNames <- c('unix', 'low', 'high', 'open', 'close', glue('{pair} volume'))
  mydata <- fromJSON(url)
  df <- as.data.frame(mydata)
  colnames(df) <- columnNames  # rename the columns
  write.csv(df, file = filename)
}


newPair <- "BTC-USD"
fileName <- glue("dailyData{newPair}.csv")
runFunc <- retreive_daily_data(newPair, filename = fileName)
runFunc



