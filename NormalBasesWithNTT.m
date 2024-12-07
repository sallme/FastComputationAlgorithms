load "GNBdata.m";
load "ENBdata.m";


/***********************************************************************
**** In this part we work with normal bases from elliptic curves. ******
************************************************************************/

/*** Preliminary Computations ***/
ListDivisor := function(n)
   ListDiv := [];
   k := n div 2;
    for d in [2..k] do
        if n mod d eq 0 then
            Append(~ListDiv, d);
        end if;
    end for;
    return ListDiv;
end function;

/******* Test if an elliptic normal basis with NTT exists *******/

ENBpossibility := function(q, n)

   possible := false;
   //nq=n whenever n is prime to q−1.
   if GCD(n, q-1) eq 1 then 
      nq := n; 
   else
      //Computation of nq using valuation
      nq := 1;
      Fn := Factorization(n);
      for L in Fn do
          l := L[1];
          if GCD(l, (q-1)) eq 1 then  
              v := L[2];
              nq := nq*(l^v);
          end if;
          if (q-1) mod l eq 0 and n mod l eq 0 then
              v := Max(2*Valuation(q-1, l)+1, 2*Valuation(n, l));
              nq := nq*(l^v);
          end if;
       end for;
   end if;
   //Main condition of existence
   if nq le Sqrt(q) then
     possible := true;
   end if;

   return possible;

end function;


/**** Works with binary extensions ****/

ENBembedDegrees := function(n, me)
   L := ListDivisor(n);
   ListDeg := [];
   for e in L do
       if e gt me then break; end if;
       if ENBpossibility(2^e, n div e) eq true then 
           Append(~ListDeg, e);      
       end if;
   end for;
   return ListDeg;
end function;



/*** Fing good elliptic embedding degree for a given range ***/

ENBrange := function(lb, ub, me)

   gd := 0;
   LowENB := [];
   for n in [lb..ub] do
       L := ENBembedDegrees(n, me);
       if #L ne 0 then 
          e := L[1];
          if e le me then 
              //printf "The order %o has a good sub degree: %o.\n", n, e;
              Append(~LowENB, <n, e>);
              gd := gd+1;
              //print L;
          end if;
       end if;
   end for;
   printf "We have %o elliptic normal bases with NTT.\n", gd;
   printf "The small embedding degrees are less than %o.\n", me;   
   printf "!!!!!!!!!!! End of computation !!!!!!!!!!!!!!.n\";
   
   return LowENB;

end function;



/**********************************************************************
** Test if a normal basis with NTT from multiplicative group exists. **
***********************************************************************/

MultNBrange := function(lb, ub, me)

    if lb mod 2 eq 0 then lb:=lb-1; end if;

    for n in [lb..ub] do
       if n mod 2 eq 0 then continue; end if;
       L := ListDivisor(n);
       if #L eq 0 or L[1] ge me+1 then continue; end if;
       Lemb := [];
          for ex in L do
              d := n div ex;        
              //The d-th roots of unity must form a proper subgroup.
              if d ge 2^ex-1 then continue; end if;
              // The embedding field must have a primitive d-th root
              if ex le me and (2^ex-1) mod d eq 0 then
                   Append(~Lemb, ex);
              end if;   
          end for;
        if #Lemb eq 0 then continue; end if;
        // Now all elements in Lemb are good embeding degrees e.
        for e in Lemb do
              //Let us find an irreducible polynomial x^d-a. 
              F2e := GF(2^e);
              a := PrimitiveElement(F2e);
        end for;
       printf "The field GF(2^%o) has a good embedding degree: %o\n", n, Lemb[1];
    end for;
   
    return "!!!!!!!!!!!!!!!!! End of Computations !!!!!!!!!!!!!!!!!";

end function;


/*********************************************************************************
*********** In this part we work with extension of normal elements ***************
*********************************************************************************/

/** Find extensions without appropriate GNBs and then work with extended bases **/
ExtNBrange := function(lb, ub, mt)
    
    /** Extensions without GNB of type<11. **/
    ExtDeg := [];
    for n in [lb..ub] do
       if n mod 8 eq 0 then 
          Append(~ExtDeg, n); 
          continue;
       end if;
       isLow := false;
       for i in [1..mt] do
           if <n, i> in LowGNB10 then 
               isLow := true; 
               break;
           end if;
       end for;
       if isLow eq false then Append(~ExtDeg, n); end if;
    end for;

    /** Now we look for Artin-Schreier extended bases **/
    ArtDeg := [];
    for n in ExtDeg do 
        if n mod 2 ne 0 then continue; end if;
        d := n div 2;
        for i in [1..mt] do
           if <d, i> in LowGNB10 then Append(~ArtDeg, <n, d, i>); end if;
        end for;
    end for;

    /** Now we look for Kummer extended bases **/
    KumDeg := [];
    for n in ExtDeg do
        if n mod 3 ne 0 then continue; end if;
        d := n div 3;
        for i in [1..mt] do
           if <d, i> in LowGNB10 then Append(~KumDeg, <n, d, i>); end if;
        end for;
    end for;
    
    /** Now we look for Witt extended bases **/
    WitDeg := [];
    for n in ExtDeg do
        if n mod 4 ne 0 then continue; end if;
        d := n div 4;
           for t in [1..mt] do
              if <d, t> in LowGNB10 then 
                  Append(~WitDeg, <n, d, t>); 
              end if;
           end for;
    end for;

    /** Let us test if an elliptic normal basis exists **/
    ElDeg := [];
    for n in ExtDeg do
        Lemb := ENBembedDegrees(n, 20);
        if #Lemb ne 0 then
            e := Lemb[1];
            if e le 20 then Append(~ElDeg, <n, e>); end if;
        end if;
    end for;

    return ArtDeg, WitDeg, KumDeg, ElDeg;   

end function;

