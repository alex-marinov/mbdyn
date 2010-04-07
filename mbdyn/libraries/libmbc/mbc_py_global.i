
%typemap(varout) unsigned *mbc_r_k_label
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_k_label_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_UINT,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_k_label
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_x
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_x_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_x
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_theta
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_theta_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_theta
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_r
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_r_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_r
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_euler_123
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_euler_123_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_euler_123
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_xp
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_xp_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_xp
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_omega
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_omega_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_omega
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_xpp
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_xpp_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_xpp
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_omegap
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_omegap_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_omegap
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) unsigned *mbc_r_d_label
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_d_label_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_UINT,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_d_label
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_f
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_f_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_f
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_r_m
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_r_m_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_r_m
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) unsigned *mbc_n_k_labels
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] =  mbc_n_k_labels_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_UINT,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_k_labels
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_x
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_x_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_x
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_theta
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_theta_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_theta
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_r
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_r_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_r
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_euler_123
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_euler_123_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_euler_123
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_xp
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_xp_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_xp
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_omega
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_omega_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_omega
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_xpp
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_xpp_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_xpp
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_omegap
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_omegap_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_omegap
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) unsigned *mbc_n_d_labels
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] =  mbc_n_d_labels_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_UINT,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_d_labels
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

%typemap(varout) double *mbc_n_f
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_f_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_f
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}


%typemap(varout) double *mbc_n_m
{
        if ($1 == NULL) {
                $result = Py_None;
        } else {
                PyArrayObject *tmp;
                int dims[1];
                dims[0] = mbc_n_m_size;
                tmp = (PyArrayObject *)PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,(char *)$1);
                $result = (PyObject *)tmp;
        }
}

%typemap(varin) double *mbc_n_m
{
        Py_INCREF($input);
        $1 = ($1_basetype *)(((PyArrayObject *)$input)->data);
}

