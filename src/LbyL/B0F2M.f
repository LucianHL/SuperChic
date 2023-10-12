      COMPLEX*16 FUNCTION B0F2M(Q2,M12)
*     ---------------------------------------
* general B0F with M12;
*
      implicit double precision(a-y)       
      COMPLEX*16 M12,LOG,SQRT,BETA
*
      BETA=SQRT(1d0-4d0*M12/Q2)
*
* B0F(-Q2;M12,M12,M12)
      B0F2M=-BETA*LOG((BETA+1d0)*(BETA+1d0)/(-4D0*M12/Q2))
*
      RETURN    
      END

    
