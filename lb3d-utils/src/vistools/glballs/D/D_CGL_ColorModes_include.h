
DISPLAYMODE_CASE_FUNCTION(D_Header_t *pHeader,D_Data_t *pArray,int colormode)
{
 switch(colormode){
  DISPLAYMODE_CASE(COLORMODE_VELOCITY,Velocity)
  DISPLAYMODE_CASE(COLORMODE_VELOCITYFIXED,VelocityFixed)
  DISPLAYMODE_CASE(COLORMODE_ANGULARVELOCITY,AngularVelocity)
  DISPLAYMODE_CASE(COLORMODE_VELOCITYNORM,VelocityNorm)
  DISPLAYMODE_CASE(COLORMODE_VELOCITYNORMFIXED,VelocityNormFixed)
  DISPLAYMODE_CASE(COLORMODE_TYPE,Type)
  DISPLAYMODE_CASE(COLORMODE_ID,Id)
  DISPLAYMODE_CASE(COLORMODE_COLOR,Color)
  DISPLAYMODE_CASE(COLORMODE_FOREGROUND,Foreground)
 }
}

//****************************************************************
// Geschwindigkeitsverteilung (Translation)

DISPLAYMODE_FUNCTION_BEGIN(Velocity)

 double VelocityMax,VelocityScale;
 double ColorValue;

 // Bestimmt den maximalen Geschwindigkeitsbetrag
 {
  double v2,v2max=T_R3_MIN;
  for(long int i=0;i<pHeader->Number;++i){
   v2=(pArray[i].vx*pArray[i].vx+pArray[i].vy*pArray[i].vy+pArray[i].vz*pArray[i].vz);
   if(v2>v2max)v2max=v2; 
  }
  VelocityMax=std::sqrt(v2max);
 }

 if(VelocityMax>0.)VelocityScale=1./VelocityMax;
 else VelocityScale=0.;

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte aus der Geschwindigkeit
 ColorValue=VelocityScale*
  std::sqrt(pArray[i].vx*pArray[i].vx+
   pArray[i].vy*pArray[i].vy+
   pArray[i].vz*pArray[i].vz);
 // Farbverlauf von Rot nach Blau
// glColor3f(ColorValue,0.,1.-ColorValue);
 glColor4f(ColorValue,0.,1.-ColorValue,0.5);

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Geschwindigkeitsverteilung (Translation)

DISPLAYMODE_FUNCTION_BEGIN(VelocityFixed)

 double VelocityScale;
 double ColorValue;

 if(m_VelocityFixed>0.)VelocityScale=1./m_VelocityFixed;
 else VelocityScale=0.;

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte aus der Geschwindigkeit
 ColorValue=VelocityScale*
  std::sqrt(pArray[i].vx*pArray[i].vx+
   pArray[i].vy*pArray[i].vy+
   pArray[i].vz*pArray[i].vz);
 if(ColorValue<1.){
  // Farbverlauf von Rot nach Blau
  glColor3f(ColorValue,0.,1.-ColorValue);
 } else {
  glColor3f(0.5,0.5,0.5);
 }

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Geschwindigkeitsverteilung (Rotation)

DISPLAYMODE_FUNCTION_BEGIN(AngularVelocity)

 double VelocityMax,VelocityScale;
 double ColorValue;

 // Bestimmt den maximalen Geschwindigkeitsbetrag
 {
  double v2,v2max=T_R3_MIN;
  for(long int i=0;i<pHeader->Number;++i){
   v2=(pArray[i].wx*pArray[i].wx+pArray[i].wy*pArray[i].wy+pArray[i].wz*pArray[i].wz);
   if(v2>v2max)v2max=v2; 
  }
  VelocityMax=std::sqrt(v2max);
 }

 if(VelocityMax>0.)VelocityScale=1./VelocityMax;
 else VelocityScale=0.;

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte aus der Geschwindigkeit
 ColorValue=VelocityScale*
  std::sqrt(pArray[i].wx*pArray[i].wx+
   pArray[i].wy*pArray[i].wy+
   pArray[i].wz*pArray[i].wz);
 // Farbverlauf von Rot nach Blau
 glColor3f(ColorValue,0.,1.-ColorValue);

DISPLAYMODE_FUNCTION_END


//****************************************************************
// Geschwindigkeitsverteilung (Translation mit Vorzugsrichtung)

DISPLAYMODE_FUNCTION_BEGIN(VelocityNorm)

 double VelocityMax,VelocityMin,VelocityScale;
 double ColorValue;

 double v;
 VelocityMin=DBL_MAX;
 VelocityMax=DBL_MIN;
 for(long int i=0;i<pHeader->Number;++i){
  v=(m_NormX*pArray[i].vx+m_NormY*pArray[i].vy+m_NormZ*pArray[i].vz);
  if(v>VelocityMax)VelocityMax=v;
  if(v<VelocityMin)VelocityMin=v;
 }
 if(fabs(VelocityMin) > fabs(VelocityMax)){
  VelocityScale=2.*fabs(VelocityMin);
 } else  {
  VelocityScale=2.*fabs(VelocityMax);
 }

 if(VelocityScale>0.)VelocityScale=1./VelocityScale;
 else VelocityScale=0.0;

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte aus der Geschwindigkeit
 ColorValue=VelocityScale*(m_NormX*pArray[i].vx+m_NormY*pArray[i].vy+m_NormZ*pArray[i].vz);
 // Farbverlauf von Rot nach Blau
 glColor3f(0.5+ColorValue,0.,0.5-ColorValue);

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Geschwindigkeitsverteilung (Translation mit Vorzugsrichtung)

DISPLAYMODE_FUNCTION_BEGIN(VelocityNormFixed)

 double VelocityScale;
 double ColorValue;

 VelocityScale=m_NormFixedScale;
 if(VelocityScale>0.)VelocityScale=1./VelocityScale;
 else VelocityScale=0.0;

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte aus der Geschwindigkeit
 ColorValue=VelocityScale*(m_NormX*pArray[i].vx+m_NormY*pArray[i].vy+m_NormZ*pArray[i].vz);
 // Farbverlauf von Rot nach Blau
 glColor3f(0.5+ColorValue,0.,0.5-ColorValue);

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Einfaerbung nach Typ

DISPLAYMODE_FUNCTION_BEGIN(Type)

DISPLAYMODE_FUNCTION_MIDDLE

 // bestimmt Farbwerte aus dem Teilchentyp
 if(m_bTypeRed && pArray[i].type==m_TypeRed)glColor3f(1.,0.,0.);
 else {
 if(m_bTypeGreen && pArray[i].type==m_TypeGreen)glColor3f(0.,1.,0.);
 else {
 if(m_bTypeBlue && pArray[i].type==m_TypeBlue)glColor3f(0.,0.,1.);
 else {
  glColor3f(0.,0.,0.);
 }}}

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Einfaerbung nach Id

DISPLAYMODE_FUNCTION_BEGIN(Id)

DISPLAYMODE_FUNCTION_MIDDLE

 // bestimmt Farbwerte aus dem Teilchentyp
 if(m_bIdRed && pArray[i].id==m_IdRed)glColor3f(1.,0.,0.);
 else {
 if(m_bIdGreen && pArray[i].id==m_IdGreen)glColor3f(0.,1.,0.);
 else {
 if(m_bIdBlue && pArray[i].id==m_IdBlue)glColor3f(0.,0.,1.);
 else {
  glColor3f(0.,0.,0.);
 }}}

DISPLAYMODE_FUNCTION_END




//****************************************************************
// Einfaerbung nach Farbwert

DISPLAYMODE_FUNCTION_BEGIN(Color)

DISPLAYMODE_FUNCTION_MIDDLE

 // berechnet Farbwerte
 glColor3f(
  0.00390625*(pArray[i].color & 255),
  0.00390625*((pArray[i].color >> 8)& 255),
  0.00390625*((pArray[i].color >> 16)& 255)
 );
DISPLAYMODE_FUNCTION_END

//****************************************************************
// Einfaerbung nach Vordergrundfarbwert

DISPLAYMODE_FUNCTION_BEGIN(Foreground)

DISPLAYMODE_FUNCTION_MIDDLE

 glColor4fv(m_paForeground);

DISPLAYMODE_FUNCTION_END












