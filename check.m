function is_OK = check( init_hio_factor,m_series,iternum,jmax )
%CHECK
assert(length(m_series)==jmax || length(iternum)~=jmax || init_hio_factor~=jmax,'Length of m_series or iternum or init_hio_factor is bad.');
for i=1:length(m_series)
    temp = m_series(i);
    assert(mod(temp,2)==0,'Numbers in m_series are not even.');
end
is_OK = 1;
end

