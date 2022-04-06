######################################################################################################################################

# All models

allinone_model <- function(data, delta, OR, endpoint, alpha=0.025){
  
  if(endpoint=="continuous"){ 
    
    m_a1 <- linear_model_a1(data, delta, alpha)
    m_a2 <- linear_model_a2(data, delta, alpha)
    m_a1_int <- linear_model_a1_int(data, delta, alpha)
    m_a2_int <- linear_model_a2_int(data, delta, alpha)
    m_b1 <- linear_model_b1(data, delta, alpha)
    m_b2 <- linear_model_b2(data, delta, alpha) 
    m_tsep <- t_test_sep(data, delta, alpha=alpha)
    m_tpol <- t_test_pol(data, delta, alpha=alpha)
    m_zsep <- z_test_sep(data, mu0=0, sigma=1, delta, alpha)
    m_zpol <- z_test_pol(data, mu0=0, sigma=1, delta, alpha)
    
    return(list(model=c("M_a1","M_a2","M_a1_int","M_a2_int","M_b1","M_b2","t_sep","t_pol","z_sep","z_pol"), 
                p_value=c(m_a1$p_val, m_a2$p_val,
                          m_a1_int$p_val, m_a2_int$p_val,
                          m_b1$p_val, m_b2$p_val,
                          m_tsep$t_test$p.value, m_tpol$t_test$p.value,
                          m_zsep$p_val, m_zpol$p_val),
                reject_h0=c(m_a1$reject_h02, m_a2$reject_h02,
                            m_a1_int$reject_h02, m_a2_int$reject_h02,
                            m_b1$reject_h02, m_b2$reject_h02,
                            m_tsep$reject_h02, m_tpol$reject_h02,
                            m_zsep$reject_h02, m_zpol$reject_h02),
                bias = c(m_a1$bias, m_a2$bias,
                         m_a1_int$bias, m_a2_int$bias, 
                         m_b1$bias, m_b2$bias, 
                         m_tsep$bias, m_tpol$bias, 
                         m_zsep$bias, m_zpol$bias)))
  }
  
  if(endpoint=="binary"){
    m_a1 <- logistic_model_a1(data, OR, alpha)
    m_a2 <- logistic_model_a2(data, OR, alpha)
    m_a1_int <- logistic_model_a1_int(data, OR, alpha)
    m_a2_int <- logistic_model_a2_int(data, OR, alpha)
    m_b1 <- logistic_model_b1(data, OR, alpha)
    m_b2 <- logistic_model_b2(data, OR, alpha)
    m_zsep <- z_prop_sep(data, OR, alpha)
    m_zpol <- z_prop_pol(data, OR, alpha)
    m_log_sep <- logistic_model_sep(data, OR, alpha)
    m_log_pol <- logistic_model_pol(data, OR, alpha)
    
    return(list(model=c("M_a1","M_a2","M_a1_int","M_a2_int","M_b1","M_b2","z_sep","z_pol", "M_log_sep", "M_log_pol"), 
                p_value=c(m_a1$p_val, m_a2$p_val,
                          m_a1_int$p_val, m_a2_int$p_val,
                          m_b1$p_val, m_b2$p_val,
                          m_zsep$p_val, m_zpol$p_val,
                          m_log_sep$p_val, m_log_pol$p_val),
                reject_h0=c(m_a1$reject_h02, m_a2$reject_h02,
                            m_a1_int$reject_h02, m_a2_int$reject_h02,
                            m_b1$reject_h02, m_b2$reject_h02,
                            m_zsep$reject_h02, m_zpol$reject_h02,
                            m_log_sep$reject_h02, m_log_pol$reject_h02),
                bias = c(m_a1$bias, m_a2$bias,
                         m_a1_int$bias, m_a2_int$bias, 
                         m_b1$bias, m_b2$bias, 
                         m_zsep$bias, m_zpol$bias,
                         m_log_sep$bias, m_log_pol$bias)))
  }
  
}

######################################################################################################################################
