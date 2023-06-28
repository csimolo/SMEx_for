INTEGER*2 FUNCTION compar(a,b)                                                                                                     
  REAL a, b                                                                                                                        
  if(a.lt.b)then                                                                                                                      
     compar = -1                                                                                                                     
  elseif(a.gt.b)then                                                                                                                  
     compar = +1                                                                                                                     
  else                                                                                                                                
     compar = 0                                                                                                                      
  endif
END FUNCTION compar
