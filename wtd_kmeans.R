library(ggplot2)
library(dplyr)
library(ggrepel)

normalize_points = function(points_clust) {
  #use feature scalingto bring range of points in the range [0,1]
  x_max = max(points_clust$x_coord);
  x_min = min(points_clust$x_coord);
  x_rng = x_max - x_min;
  
  y_max = max(points_clust$y_coord);
  y_min = min(points_clust$y_coord);
  y_rng = y_max - y_min;
  
  points_clust = mutate(points_clust, x_coord_n = (x_coord - x_min) / x_rng, y_coord_n = (y_coord - y_min) / y_rng);
  
  return(points_clust);
}

calc_dist = function(x1, y1, x2, y2) {
  return(sqrt((x1-x2)^2 + (y1-y2)^2));
}

calc_centroids = function(points_clust) {
  #j = cluster j
  #i = point i in cluster j
  #N = number of customers in cluster j
  #w_i = weight for point i
  #wtd_centroid_j = weight x,y coords of centroid for cluster j
  #               = \sum_{i=1}^{N} (x_i,y_i) * w_i / \sum{i=1}^{N} w_i
  #wtd mse = E(J) = Sum ( all X(I) in cluster J ) W(I) * || X(I) - Z(J) ||^2 / \sum{i=1}^{N} w_i
  d = summarise(group_by(points_clust, cluster),
                wtd_x_coord = sum(x_coord * weight),
                wtd_y_coord = sum(y_coord * weight),
                wtd_x_coord_n = sum(x_coord_n * weight),
                wtd_y_coord_n = sum(y_coord_n * weight),
                tot_wtds = sum(weight), na.rm = TRUE);
  d = mutate(d,
             c_x_coord = wtd_x_coord / tot_wtds, 
             c_y_coord = wtd_y_coord / tot_wtds,
             c_x_coord_n = wtd_x_coord_n / tot_wtds, 
             c_y_coord_n = wtd_y_coord_n / tot_wtds);
  
  d = mutate(d, cluster = as.integer(cluster));
  points_clust = mutate(points_clust, cluster = as.integer(cluster));
  
  #calc wtd mse for each cluster with centroid for cluster
  e = left_join(points_clust, d, by = c("cluster" = "cluster")) %>%
      mutate(dist = calc_dist(x_coord_n, y_coord_n, c_x_coord_n, c_y_coord_n));
  
  e = summarise(group_by(e, cluster), sse = sum(dist * weight) / sum(weight));

  
  centroids = left_join(d, e, by = c("cluster" = "cluster")) %>% 
              select(cluster, c_x_coord, c_y_coord, c_x_coord_n, c_y_coord_n, sse);
  
  return(centroids);
}

assign_points = function(points_clust, centroids) {
  #assignment based on normalized distance
  d = merge(points_clust, centroids, by = NULL) %>%
      mutate(dist_clust = calc_dist(x_coord_n, y_coord_n, c_x_coord_n, c_y_coord_n)) %>%
      group_by(prem_id_org) %>% slice(which.min(dist_clust));
  
  d = mutate(d, chg = !(cluster.x == cluster.y)) %>%
      select(prem_id_org, x_coord, y_coord, weight, cluster = as.integer(cluster.y), 
             x_coord_n, y_coord_n, chg);
  
  #d = rename(d, cluster = cluster.y);

  return(d)
}

plot_points = function(points_clust, centroids, i, scale = "o") {
  
  if(scale == "o") {
    plot = ggplot() +
           geom_point(data = points_clust, aes(x=x_coord, y=y_coord, col = factor(cluster)), alpha = 0.2) +
           geom_point(data = centroids, aes(x=c_x_coord, y = c_y_coord), shape = 4, size = 3, col = "black") +  
           coord_cartesian(xlim=c(-8000,28000), ylim=c(0, 45)) +
           labs(title = paste("Iteration",i));
  }
  
  if(scale == "n"){
    plot = ggplot() +
      geom_point(data = points_clust, aes(x=x_coord_n, y=y_coord_n, col = factor(cluster)), alpha = 0.2) +
      geom_point(data = centroids, aes(x=c_x_coord_n, y = c_y_coord_n), shape = 4, size = 3, col = "black") +  
      coord_cartesian(xlim=c(-0.1,1.1), ylim=c(-0.1, 1.1)) +
      labs(title = paste("Original Scale - Iteration",i));
    
    plot = plot + geom_text_repel(aes(label=sse, x=c_x_coord_n, y=c_y_coord_n), size=4, data=centroids)     + 
      theme(legend.position = "None")   # text;
    
  }
  print(plot);
}

main = function(points, k, max_itr) {
  points_clust = mutate(points, cluster = as.integer(NA))
  
  #normalize points
  points_clust = normalize_points(points_clust);
  
  #select initial starting centroids at random
  set.seed(Sys.time(), kind = NULL, normal.kind = NULL)
  x_init = runif(k);
  y_init = runif(k);
  centroids = data.frame(cluster = as.integer(1:k), c_x_coord = NA, c_y_coord = NA,
                         c_x_coord_n = x_init, c_y_coord_n = y_init);

  
  #initial assignment
  i = 1;
  points_clust = assign_points(points_clust, centroids);
  centroids = calc_centroids(points_clust);
  print(paste("Iteration ", i, ": ", c(centroids$sse, sum(centroids$sse)), sep=""));
  plot_points(points_clust, centroids, i, "n");
  
  done = FALSE;
  while (done == FALSE && i < max_itr) {
    i = i + 1;
    points_clust = assign_points(points_clust, centroids);
    if (max(points_clust$chg) == FALSE) {
      done = TRUE
    } else {
      centroids = calc_centroids(points_clust);
      #print(paste("Iteration ", i, ": ", c(centroids$sse, sum(centroids$sse)), sep=""));
      #plot_points(points_clust, centroids, i, "n");
    }
  }
  print(paste("Iteration ", i, ": ", c(centroids$sse, sum(centroids$sse)), sep=""));
  plot_points(points_clust, centroids, i, "n");
  plot_points(points_clust, centroids, i, "o");
  
  #K-means alg:
  #select k points as the initial centroids
  #repeat
  #  form k clusters by assigning each point to its closest centroid
  #  recompute the centroid of each cluster
  #until centroids do not change
  
  #Inputs
  #k - integer: number of clusters
  #max_itr - integer: maximum number of iterations
  #points - dataframe: prem_id_org, x_coord, y_coord, weight
  #Return
  #points_clust - dataframe: prem_id_org, x_coord, y_coord, weight, cluster, x_coord_n, y_coord_n
  #centroids - dataframe: cluster, c_x_coord, c_y_coord, c_x_coord_n, c_y_coord_n, sse
  
  return(list(points_clust, centroids));
}

k = 4;
max_itr = 50;
points = rename(df_nem, x_coord = ann_net_usg, y_coord = ann_pk_del_dmd, weight = wtd);

print(paste("K =",k));
start_time = Sys.time();
out = main(points, k, max_itr);
end_time = Sys.time();
print(end_time - start_time);

