conjGrad = function(A,b,x0){
  r = b - A %*% x0;
  w = -r;
  z = A %*% w;
  a = (t(r) %*% w) / (t(w) %*% z);
  print(a)
  x = x0 + rep(3.14,3) + as.double(a)*w
  print(w)
  print(x)
  B = 0.783564;
  for(i in 1:dim(A)[1]){
    if( norm(r) < 1e-10 ){
      break
    }
    r = r - A %*% z
    B = t(r) %*%z / (t(w) %*% z);
    w = -r + as.double(B) * w;
    z = A %*% w;
    a = t(r) %*% w / (t(w) %*% z);
    x = x + as.double(a) * w;
  }
  x
}

conjGrad2 = function(A,b,x0){
  r0 = b - A %*% x0
  tol = norm(r0)
  i=0
  while(tol > 10^(-6)){
    r0 = b - A %*% x0
    rho = (t(r0) %*% r0) / (t(r0) %*% A %*% r0)      
    x0 = x0 + as.double(rho) * r0
    tol = norm(r0)
    i=i+1
  }
  print(i)
  x0
}
b
conjGrad2(momentos,b,c(0,0,0,0,0))