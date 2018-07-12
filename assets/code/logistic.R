#########################################################################
# A simple plotting procedure for May's logistic equation of population #
# growth.                                                               #
#                                                                       #
# by Mario dos Reis. -- October 2003                                    #
#########################################################################

#########################################################################
# Chaos in a simple population system
#########################################################################

# non linear, chaotic function
# r is growth rate and x is population size (from 0 to 1)
lg = function(x, r) r * x * (1 - x)


gen = 50       # number of initial generations
r.min = 3.2
r.max = 4
p = '.'
cl = 'black'   # plotting color
div = 4000     # plot resolution along the x axis
k = 16         # generations to be plotted, increase to 64 to get a high
               # resolution plot

# full plot
plot(x=0, y=0, type = 'n', xlim = c(r.min, r.max), ylim = c(0, 1),
     xlab = 'r', ylab = 'X')

for(r in seq(r.min, r.max, by = r.max/div)) {

  x.init = .01
  
  for(i in seq(0, gen)) {
   
    x.next = lg(x.init, r)
    x.init = x.next
    
   }

  for(i in seq(0, k)) {

    x.next = lg(x.init, r)
    x.init = x.next
    
    points(r, x.next, pch = p, col = cl)

  }
}

