function is_OK = check( m_series,iternum,jmax )
%CHECK
assert(length(m_series)==jmax || length(iternum)~=jmax,'Length of m_series is bad.');
for i=1:length(m_series)
    temp = m_series(i);
    assert(mod(temp,2)==0,'Numbers in m_series are not even.');
end
is_OK = 1;
end

