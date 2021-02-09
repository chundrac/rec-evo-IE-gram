require(rstan)
require(phytools)

#feature indices {1..65}
feat.ind = as.integer(commandArgs(trailingOnly=TRUE)[1])
print(feat.ind)

if (!dir.exists('stanfits/')) {
	dir.create('stanfits/')
}

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)

IE.trees <- read.newick('IE.newick')

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

model.code <- "data {
  int<lower=1> N; //number of tips+internal nodes+root
  int<lower=1> T; //number of tips
  int<lower=1> B; //number of branches
  int<lower=1> F;                       //number of states
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];                //length of each branch
  int<lower=0,upper=1> tiplik[T,F];     //likelihoods for data at tips in tree
  }
parameters {
  real<lower=0> R[F*(F-1)];             //rates to be put into matrix
  }
transformed parameters {
  matrix[F,F] Q;                       //rate matrix
  //fill off-diagonal cells of rate matrix
  {
    int k;
    k = 1;
    for (i in 1:F) {
      for (j in 1:F) {
        if (i != j) {
          Q[i,j] = R[k];
          k = k + 1;
        }
      }
    }
  }
  //fill diagonal rows of rate matrix
  for (i in 1:F) {
    real z;
    z = 0;
    for (j in 1:F) {
      if (i != j) {
        z = z - Q[i,j];
      }
    }
    Q[i,i] = z;
  }
}
model {
  matrix[N,F] lambda; //matrix of likelihoods for root, internal nodes, and tips
  matrix[F,F] pi_matrix;                //holds approximation of stationary probabilities in each row
  vector[F] pi;                         //stationary probability
  for (i in 1:F*(F-1)) {
    R[i] ~ gamma(1,1);                //shape/rate parameterization
  }
  for (t in 1:T) {
    for (f in 1:F) {
      lambda[t,f] = log(tiplik[t,f]);
    }
  }
  for (n in (T+1):N) {
    for (f in 1:F) {
      lambda[n,f] = 0;
    }
  }
  for (b in 1:B) {
    matrix[F,F] P;
    P = matrix_exp(brlen[b]*Q); //via matrix exponentiation
    for (f in 1:F) {
      lambda[parent[b],f] = lambda[parent[b],f] + log(dot_product(P[f],exp(lambda[child[b]])));
    }
  }
  pi_matrix = matrix_exp(1000*Q);     //
  for (f in 1:F) {                      //dirty way of approximating
    pi[f] = pi_matrix[1,f];             //stationary probabilities of
  }                                     //each character state;
  target += log(dot_product(pi,exp(lambda[parent[B]])));   //increment log posterior by tree likelihood
}
generated quantities {
  matrix[N,F] lambda;
  matrix[F,F] pi_matrix;                //holds approximation of stationary probabilities in each row
  vector[F] pi;                         //stationary probability
  vector[F] psi;
  vector[F] phi;
  int z[F];
  for (t in 1:T) {
    for (f in 1:F) {
      lambda[t,f] = log(tiplik[t,f]);
    }
  }
  for (n in T+1:N) {
    for (f in 1:F) {
      lambda[n,f] = 0;
    }
  }
  for (b in 1:B) {
    matrix[F,F] P;
    P = matrix_exp(brlen[b]*Q); //via matrix exponentiation
    for (f in 1:F) {
      lambda[parent[b],f] = lambda[parent[b],f] + log(dot_product(P[f],exp(lambda[child[b]])));
    }
  }
  pi_matrix = matrix_exp(1000*Q);       //
  for (f in 1:F) {                      //dirty way of approximating
    pi[f] = pi_matrix[1,f];             //stationary probabilities of
  }                                     //each character state;
  for (f in 1:F) {
    psi[f] = lambda[parent[B],f]+log(pi[f]);
    }
  for (f in 1:F) {
    phi[f] = exp(psi[f]-log_sum_exp(psi));
    }
  z = multinomial_rng(phi,1);
}"

print(levels(diacl.data[,feat.ind]))

likelihood <- to.matrix(as.factor(diacl.data[,feat.ind]),seq=levels(as.factor(diacl.data[,feat.ind])))
rownames(likelihood) <- rownames(diacl.data)
likelihood[rowSums(likelihood)==0,] <- likelihood[rowSums(likelihood)==0,]+1
fit.list <- list()	
for (i in 1:length(IE.trees)) {
    tree <- IE.trees[[i]]
   	tree <- reorder.phylo(tree,'pruningwise')
    curr.states <- likelihood[tree$tip.label,]
   	parent <- tree$edge[,1]
   	child <- tree$edge[,2]
   	b.lens <- tree$edge.length
   	N <- length(unique(c(parent,child)))
   	T <- length(child[which(!child %in% parent)])
       data.list <- list(N=N,
             T=T,
             B=length(parent),
             brlen=b.lens/1000,
             child=child,
             parent=parent,
             tiplik=curr.states,
             F=ncol(curr.states))
       fit <- stan(model_code=model.code,data=data.list,chains=3,thin=10)
   	fit.list[[i]] <- fit
	fit.full <- sflist2stanfit(fit.list)
	save.image(file=paste('stanfits/',colnames(diacl.data)[feat.ind],'.Rdata',sep=''))
   }


