
#[allow(dead_code)]
pub mod kali_vector {

    use std::ops;

    #[derive(Debug)]
    pub enum DivisionError {
        DivisionByZero
    }

    #[derive(Debug)]
    pub struct KaliVector {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    impl KaliVector {

        pub fn new(x: f64, y: f64, z: f64) -> KaliVector {
            KaliVector { x: x, y: y, z: z }
        }

        pub fn new_zero() -> KaliVector {
            KaliVector { x: 0.0, y: 0.0, z: 0.0 }
        }

        pub fn dot_product(self, vec: KaliVector) -> f64 {
            self.x * vec.x +
            self.y * vec.y +
            self.z * vec.z
        }

        pub fn cross_product(self, vec: KaliVector) -> KaliVector {
            KaliVector {
                x: self.y * vec.z - self.z * vec.y,
                y: self.z * vec.x - self.x * vec.z,
                z: self.x * vec.y - self.y * vec.x
            }
        }

        pub fn magnitude(self) -> f64 {
            let radicand: f64 = self.x * self.x + 
                self.y * self.y + 
                self.z * self.z;
            radicand.sqrt()
        }

        pub fn normalize(&mut self) {
            let vec = self.clone();
            let magnitude = vec.magnitude();

            if magnitude > 0.0 {
                let one_over_mag = 1.0 / magnitude;

                *self = Self {
                    x: vec.x * one_over_mag,
                    y: vec.y * one_over_mag,
                    z: vec.z * one_over_mag
                };
            }
        }
    }

    impl Copy for KaliVector { }

    impl Clone for KaliVector {
        fn clone(&self) -> KaliVector {
            *self
        }
    }

    impl ops::Add<KaliVector> for KaliVector {
        type Output = Self;

        fn add(self, _rhs: Self) -> Self {
            Self {
                x: self.x + _rhs.x,
                y: self.y + _rhs.y,
                z: self.z + _rhs.z
            }
        }
    }

    impl ops::AddAssign for KaliVector {
        fn add_assign(&mut self, _rhs: Self) {
            *self = Self {
                x: self.x + _rhs.x,
                y: self.y + _rhs.y,
                z: self.z + _rhs.z
            };
        }
    }

    impl ops::Sub<KaliVector> for KaliVector {
        type Output = Self;

        fn sub(self, _rhs: Self) -> Self {
            Self {
                x: self.x - _rhs.x,
                y: self.y - _rhs.y,
                z: self.z - _rhs.z
            }
        }
    }

    impl ops::SubAssign for KaliVector {
        fn sub_assign(&mut self, _rhs: Self) {
            *self = Self {
                x: self.x - _rhs.x,
                y: self.y - _rhs.y,
                z: self.z - _rhs.z
            };
        }
    }

    impl ops::Mul<f64> for KaliVector {
        // TODO: change to T: Float
        type Output = Self;

        fn mul(self, scalar: f64) -> Self {
            Self {
                x: self.x * scalar,
                y: self.y * scalar,
                z: self.z * scalar
            }
        }
    }

    impl ops::MulAssign<f64> for KaliVector {
        fn mul_assign(&mut self, scalar: f64) {
            *self = Self {
                x: self.x * scalar,
                y: self.y * scalar,
                z: self.z * scalar,
            };
        }
    }

    impl ops::Div<f64> for KaliVector {
        type Output = Self;

        fn div(self, scalar: f64) -> Self {
            if scalar == 0.0 {
                panic!("divide by zero");
            }
            
            KaliVector {
                x: self.x / scalar,
                y: self.y / scalar,
                z: self.z / scalar
            }
        }
    }

    impl ops::DivAssign<f64> for KaliVector {

        fn div_assign(&mut self, scalar: f64) {
            if scalar == 0.0 {
                panic!("divide by zero");
            }

            *self = Self {
                x: self.x / scalar,
                y: self.y / scalar,
                z: self.z / scalar
            };
        }
    }
}



#[allow(dead_code)]
pub mod kali_matrix {
    
    use std::ops;
    
    #[derive(Debug)]
    pub struct KaliMatrix {
        pub values: [f64; 16]
    }
    
    impl KaliMatrix {
        
        pub fn from_array(arr: [f64; 16]) -> KaliMatrix {
            KaliMatrix {
                values: arr
            }
        }
        
        pub fn new_zero() -> KaliMatrix {
            let arr: [f64; 16] = [0.0; 16];
            KaliMatrix {
                values: arr
            }
        }
        
        pub fn scale(&mut self, s: f64) {
            for i in 0..=15 {
                self.values[i] *= s;
            }
        }
    }

    impl Copy for KaliMatrix { }

    impl Clone for KaliMatrix {
        fn clone(&self) -> KaliMatrix {
            *self
        }
    }
    
    /* column major format
    0   1   2   3
    4   5   6   7
    8   9   10  11
    12  13  14  15
    */
    
    impl ops::Add for KaliMatrix {
        type Output = Self;
        
        fn add(self, _rhs: Self) -> KaliMatrix {
            let left = self.values;
            let right  = _rhs.values;
            
            let mut out = KaliMatrix::new_zero();
            
            out.values[0] = left[0] + right[0];
            out.values[4] = left[4] + right[4];
            out.values[8] = left[8] + right[8];
            out.values[12] = left[12] + right[12];
            
            out.values[1] = left[1] + right[1];
            out.values[5] = left[5] + right[5];
            out.values[9] = left[9] + right[9];
            out.values[13] = left[13] + right[13];
            
            out.values[2] = left[2] + right[2];
            out.values[6] = left[6] + right[6];
            out.values[10] = left[10] + right[10];
            out.values[14] = left[14] + right[14];
            
            out.values[3] = left[3] + right[3];
            out.values[7] = left[7] + right[7];
            out.values[11] = left[11] + right[11];
            out.values[15] = left[15] + right[15];
            
            out
        }
    }
}


#[cfg(test)]
mod kali_vector_tests {
    use crate::kali_vector::KaliVector;

    #[test]
    fn test_kali_vector_new_zero() {
        let vec = KaliVector::new_zero();
        assert_eq!(vec.x, 0.0);
        assert_eq!(vec.y, 0.0);
        assert_eq!(vec.z, 0.0);
    }

    #[test]
    fn test_kali_vector_new() {
        let vec = KaliVector::new(0.0, 1.0, 2.0);
        assert_eq!(vec.x, 0.0);
        assert_eq!(vec.y, 1.0);
        assert_eq!(vec.z, 2.0);
    }

    #[test]
    fn test_adding_kali_vectors() {
        let vec1 = KaliVector::new(1.0, 2.0, 3.0);
        let vec2 = KaliVector::new(1.0, 2.0, 3.0);
        let vec3 = vec1 + vec2;
        assert_eq!(vec3.x, 2.0);
        assert_eq!(vec3.y, 4.0);
        assert_eq!(vec3.z, 6.0);
    }

    #[test]
    fn test_add_assign_overload() {
        let mut vec1 = KaliVector::new(1.0, 2.0, 3.0);
        let vec2 = KaliVector::new(1.0, 2.0, 3.0);
        vec1 += vec2;

        assert_eq!(vec1.x, 2.0);
        assert_eq!(vec1.y, 4.0);
        assert_eq!(vec1.z, 6.0);
    }

    #[test]
    fn test_sub_two_kali_vectors() {
        let vec1 = KaliVector::new(1.0, 2.0, 3.0);
        let vec2 = KaliVector::new(1.0, 2.0, 3.0);
        let vec3 = vec1 - vec2;
        assert_eq!(vec3.x, 0.0);
        assert_eq!(vec3.y, 0.0);
        assert_eq!(vec3.z, 0.0);
    }

    #[test]
    fn test_sub_assign() {
        let mut vec1 = KaliVector::new(1.0, 2.0, 3.0);
        let vec2 = KaliVector::new(1.0, 2.0, 3.0);
        vec1 -= vec2;

        assert_eq!(vec1.x, 0.0);
        assert_eq!(vec1.y, 0.0);
        assert_eq!(vec1.z, 0.0);
    }

    #[test]
    fn test_multiply_kali_vector_by_scalar() {
        let vec1 = KaliVector::new_zero();
        let vec2 = vec1 * 1.0;
        assert_eq!(vec2.x, 0.0);
        assert_eq!(vec2.y, 0.0);
        assert_eq!(vec2.z, 0.0);
    }

    #[test]
    fn test_multiply_assign() {
        let mut vec1 = KaliVector::new(1.0, 1.0, 1.0);
        vec1 *= 3.0;
        assert_eq!(vec1.x, 3.0);
        assert_eq!(vec1.y, 3.0);
        assert_eq!(vec1.z, 3.0);
    }

    #[test]
    fn test_divide_kali_vector_by_scalar() {
        let vec1 = KaliVector::new(4.0, 4.0, 4.0);
        let vec2 = vec1 / 2.0;

        assert_eq!(vec2.x, 2.0);
    }

    #[test]
    #[should_panic]
    #[allow(unused_variables)]
    fn test_divide_by_zero_on_kali_vector() {
        let vec1 = KaliVector::new(4.0, 4.0, 4.0);
        let vec2 = vec1 / 0.0;
    }

    #[test]
    fn test_divide_assign_kali_vector() {
        let mut vec1 = KaliVector::new(4.0, 4.0, 4.0);
        vec1 /= 2.0;
        assert_eq!(vec1.x, 2.0);
    }

    #[test]
    #[should_panic]
    fn test_divide_assign_panic_on_divide_by_zero() {
        let mut vec1 = KaliVector::new(4.0, 4.0, 4.0);
        vec1 /= 0.0;
    }

    #[test]
    fn test_dot_product() {
        let vec1 = KaliVector::new(2.0, 3.0, 1.0);
        let vec2 = KaliVector::new(1.0, 2.0, 0.0);
        let dot_product = vec1.dot_product(vec2);
        assert_eq!(dot_product, 8.0);
    }

    #[test]
    fn test_cross_product() {
        let vec = KaliVector::new(2.0, 3.0, 1.0);
        let vec2 = KaliVector::new(1.0, 2.0, 0.0);

        let final_vec = vec.cross_product(vec2);
        assert_eq!(final_vec.x, -2.0);
        assert_eq!(final_vec.y, 1.0);
        assert_eq!(final_vec.z, 1.0);
    }

    #[test]
    fn test_magnitue() {
        let vec = KaliVector::new(2.0, 3.0, 1.0);
        let mut magnitude = vec.magnitude();
        magnitude = (magnitude * 10000.0).round() / 10000.0;
        assert_eq!(magnitude, 3.7417);
    }

    #[test]
    fn test_normalize() {
        let mut vec = KaliVector::new(2.0, 3.0, 1.0);

        vec.normalize();
        
        assert_ne!(vec.x, 2.0);
        assert_ne!(vec.y, 3.0);
        assert_ne!(vec.z, 1.0);
    }
}

#[cfg(test)]
mod kali_matrix_tests {
    use crate::kali_matrix::KaliMatrix;
    
    #[test]
    fn test_new_kali_matrix_from_array() {
        let arr: [f64; 16] = [1.0; 16];
        let matrix = KaliMatrix::from_array(arr);
        for i in 0..=15 {
            assert_eq!(matrix.values[i], 1.0);
        }
    }

    #[test]
    fn test_new_zeroed_kali_matrix() {
        let matrix = KaliMatrix::new_zero();
        for i in 0..=15 {
            assert_eq!(matrix.values[i], 0.0);
        }
    }

    #[test]
    fn test_kali_matrix_scale() {
        let arr: [f64; 16] = [1.0; 16];
        let mut matrix = KaliMatrix::from_array(arr);

        matrix.scale(5.0);
        println!("{:?}\n", matrix);
        for i in 0..=15 {
            assert_eq!(matrix.values[i], 5.0);
        }
    }
}
