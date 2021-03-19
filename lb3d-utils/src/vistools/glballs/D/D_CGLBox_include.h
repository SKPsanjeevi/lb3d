int D_CGLBox::getDisplayMode(void) const
{
 return m_CGL.m_DisplayMode;
}

int D_CGLBox::getColorMode(void) const
{ 
 return m_CGL.m_ColorMode;
}

//
/*
*/
int D_CGLBox::getXRotation(void) const
{
 return ((int)m_CGL.m_xRot) % 360;
}

//
/*
*/
int D_CGLBox::getYRotation(void) const
{
 return ((int)m_CGL.m_yRot) % 360;
}

//
/*
*/
int D_CGLBox::getZRotation(void) const
{
 return ((int)m_CGL.m_zRot) % 360;
}

//
/*
*/
int D_CGLBox::getScale(void) const
{
 return (int)(std::log(m_CGL.m_scale/0.2)/0.03);
}

//
/*
*/
int D_CGLBox::getViewScale(void) const
{
 return (int)(std::log(m_CGL.m_viewscale/0.2)/0.03);
}

//
/*
*/
int D_CGLBox::getOffset(void) const
{
 return (int)(10*m_CGL.m_offset);
}

//
/*
*/
//int D_CGLBox::getOffsetX(void) const
//{
//return ;
//}

//
/*
*/
//int D_CGLBox::getOffsetY(void) const
//{
//return ;
//}

