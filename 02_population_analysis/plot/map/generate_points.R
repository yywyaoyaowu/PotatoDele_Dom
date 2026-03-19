library(sf)
library(rnaturalearth)
library(dplyr)
library(units)

Args <- commandArgs(TRUE)
mapping <- Args[1]
num <- Args[2]
num <- as.numeric(num)

generate_points <- function(country, n, buffer_distance_km = 500, max_attempts = 100, seed = 123) {
  set.seed(seed)
  
  Country <- ne_countries(country = country, scale = "medium", returnclass = "sf")
  Country <- st_make_valid(Country)
  
  original_crs <- st_crs(Country)
  
  centroid <- st_centroid(Country)
  centroid_coords <- st_coordinates(centroid)
  utm_zone <- floor((centroid_coords[1] + 180) / 6) + 1
  utm_crs <- ifelse(centroid_coords[2] >= 0, 
                   paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84"),
                   paste0("+proj=utm +zone=", utm_zone, " +south +datum=WGS84"))
  
  Country_utm <- st_transform(Country, utm_crs)
  
  buffer_dist_m <- buffer_distance_km * 1000
  points_sf <- NULL
  attempt <- 1
  
  while(attempt <= max_attempts && is.null(points_sf)) {
    tryCatch({
      buffered_country <- st_buffer(Country_utm, dist = -buffer_dist_m)
      
      if (st_is_valid(buffered_country) && st_area(buffered_country) > set_units(0, "m^2")) {
        buffered_country_wgs84 <- st_transform(buffered_country, 4326)
        
        points_sf <- st_sample(buffered_country_wgs84, size = n, type = "random")
        
        if (length(points_sf) == n) {
          break
        }
      }
    }, error = function(e) {
      buffer_dist_m <- buffer_dist_m * 0.7
      attempt <- attempt + 1
    })
  }
  
  if (is.null(points_sf) || length(points_sf) < n) {
    warning("缓冲方法失败或生成点数不足，使用原始采样方法并添加边界距离检查")
    points_sf <- st_sample(Country, size = n * 2, type = "random")
    
    country_boundary <- st_boundary(Country)
    distances <- st_distance(points_sf, country_boundary)
    
    min_distance <- set_units(10, "km")
    valid_points <- which(as.numeric(distances) > as.numeric(min_distance))
    
    if (length(valid_points) >= n) {
      points_sf <- points_sf[valid_points[1:n]]
    } else {
      warning("无法找到足够多的远离边界的点，返回所有可用点")
      points_sf <- points_sf[valid_points]
    }
  }
  
  points_sf <- st_transform(points_sf, crs = 4326)
  coords <- st_coordinates(points_sf)
  colnames(coords) <- c("longitude", "latitude")
  
  return(as.data.frame(coords))
}

points_df <- generate_points(country = mapping, n = num, buffer_distance_km = 100)
colnames(points_df) <- c("Longitude", "Latitude")
points_df <- points_df[, c(2, 1)]
print(points_df, row.names = FALSE)



