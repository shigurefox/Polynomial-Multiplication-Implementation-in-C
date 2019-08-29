/*
 * poly.c
 *
 *  Created on: Aug 27, 2019
 *      Author: shigurefox
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "strman.h"
#include "poly.h"

// ========================================
//   Modular Multiplicative Inverse Table
// ========================================

const int mod_inv[MODULUS] = {0, 1, 2296, 3061, 1148, 3673, 3826, 656, 574, 4081, \
		4132, 3339, 1913, 2119, 328, 4285, 287, 4321, 4336, 2658, \
		2066, 1749, 3965, 3593, 3252, 2571, 3355, 4421, 164, 2533, \
		4438, 1481, 2439, 1113, 4456, 3804, 2168, 1489, 1329, 3767, \
		1033, 112, 3170, 1388, 4278, 4489, 4092, 2149, 1626, 3373, \
		3581, 4501, 3973, 693, 4506, 1586, 82, 886, 3562, 3346, \
		2219, 1430, 3036, 583, 3515, 1342, 2852, 3015, 2228, 2728, \
		1902, 194, 1084, 4025, 3040, 857, 2960, 477, 4179, 2034, \
		2812, 3004, 56, 3706, 1585, 4537, 694, 3905, 2139, 3972, \
		4540, 3582, 2046, 2024, 3370, 2368, 813, 142, 3982, 371, \
		4086, 4091, 4546, 4279, 4282, 1268, 2642, 3218, 2253, 2822, \
		793, 3557, 41, 1422, 443, 2555, 1781, 2786, 1673, 1929, \
		3405, 1973, 715, 3098, 1518, 4187, 2587, 723, 4053, 1993, \
		671, 3820, 1426, 3659, 3803, 4557, 1114, 1508, 1364, 3435, \
		951, 3777, 97, 610, 542, 2343, 4308, 4185, 1520, 493, \
		2724, 3010, 1480, 4561, 2534, 3969, 4385, 2778, 1017, 231, \
		1406, 1825, 1502, 169, 28, 2059, 1853, 3079, 3088, 163, \
		4564, 3356, 347, 2468, 4248, 1679, 3365, 4176, 1986, 1385, \
		2270, 2156, 1791, 2007, 1023, 1216, 1012, 2897, 1685, 3255, \
		1184, 2620, 2702, 2688, 71, 3508, 1991, 4055, 2481, 1638, \
		2043, 1005, 4341, 4297, 2273, 2777, 4435, 3970, 2141, 659, \
		634, 2611, 1321, 1595, 1609, 2114, 3422, 2835, 1411, 2872, \
		2692, 1745, 4074, 1647, 2316, 1816, 711, 1800, 2517, 2085, \
		3573, 159, 3186, 4059, 1393, 1348, 3132, 678, 3260, 826, \
		3998, 381, 3282, 4062, 2653, 2511, 1549, 3736, 759, 4296, \
		4389, 1006, 3589, 744, 2657, 4573, 4322, 4198, 3292, 3492, \
		2631, 2832, 1910, 1798, 713, 1975, 4125, 1324, 4197, 4335, \
		4574, 288, 557, 1194, 754, 3990, 682, 895, 4013, 2205, \
		2771, 3431, 4184, 4445, 2344, 3850, 305, 16, 271, 4035, \
		3467, 3108, 2154, 2272, 4388, 4342, 760, 1654, 2542, 1689, \
		1362, 1510, 1505, 2894, 740, 286, 4576, 329, 1267, 4487, \
		4280, 4281, 4488, 4547, 1389, 1953, 2804, 420, 2411, 1065, \
		703, 2603, 3208, 3127, 751, 3023, 2380, 2471, 14, 307, \
		3325, 1068, 3222, 2716, 3835, 603, 1544, 940, 2377, 474, \
		2282, 552, 1678, 4417, 2469, 2382, 1234, 172, 2124, 1105, \
		3135, 2459, 3978, 2289, 2088, 957, 993, 643, 2988, 1752, \
		1135, 4006, 1078, 2188, 3191, 805, 3299, 1964, 2807, 2563, \
		608, 99, 506, 1477, 3744, 2926, 3138, 548, 3923, 1708, \
		592, 241, 1310, 2757, 1351, 2850, 1344, 3725, 2331, 1133, \
		1754, 3182, 3291, 4334, 4323, 1325, 3536, 2070, 819, 2750, \
		3317, 2576, 2798, 2586, 4466, 1519, 4444, 4309, 3432, 449, \
		3684, 2033, 4513, 478, 1985, 4414, 3366, 1145, 2625, 1786, \
		317, 2181, 3601, 1259, 2956, 3662, 3093, 1516, 3100, 3264, \
		1057, 2418, 1711, 3817, 3713, 781, 3001, 914, 1436, 3723, \
		1346, 1395, 3168, 114, 2037, 3549, 3119, 3225, 1158, 409, \
		908, 2097, 2651, 4064, 900, 3471, 3554, 874, 3338, 4581, \
		4082, 3625, 2375, 942, 1593, 1323, 4325, 1976, 2992, 3710, \
		674, 926, 1566, 2213, 339, 2310, 1630, 77, 413, 2607, \
		1999, 1174, 2486, 3669, 1641, 2783, 2031, 3686, 3622, 3117, \
		3551, 1468, 3070, 149, 1868, 3747, 2675, 3307, 2148, 4545, \
		4490, 4087, 503, 502, 4090, 4491, 372, 3115, 3624, 4131, \
		4582, 575, 2161, 2649, 2099, 1774, 1646, 4369, 1746, 2353, \
		3611, 3701, 1416, 2256, 955, 2090, 899, 4138, 2652, 4348, \
		3283, 1392, 4358, 3187, 662, 2480, 4394, 1992, 4463, 724, \
		2287, 3980, 144, 2249, 2574, 3319, 597, 1452, 377, 669, \
		1995, 2308, 341, 2914, 2743, 3466, 4302, 272, 3398, 2932, \
		3681, 2496, 4011, 897, 2092, 3039, 4518, 1085, 1172, 2001, \
		1925, 3393, 2448, 3934, 8, 511, 2431, 2204, 4313, 896, \
		4029, 2497, 1554, 63, 1077, 4230, 1136, 962, 2194, 2011, \
		2171, 2882, 380, 4351, 827, 1304, 1271, 546, 3140, 3449, \
		681, 4316, 755, 335, 3048, 2231, 1447, 2027, 370, 4493, \
		143, 4050, 2288, 4239, 2460, 2456, 2929, 692, 4539, 4502, \
		2140, 4384, 4436, 2535, 2244, 3592, 4569, 1750, 2990, 1978, \
		3272, 2401, 1402, 3931, 210, 1981, 3501, 3791, 2828, 2062, \
		2647, 2163, 3597, 357, 1604, 2235, 3859, 1703, 2671, 2391, \
		3807, 945, 1190, 3902, 3531, 764, 7, 4018, 2449, 209, \
		3958, 1403, 534, 2112, 1611, 1650, 1358, 1707, 4213, 549, \
		2597, 130, 772, 880, 470, 3666, 3484, 1458, 237, 1332, \
		1141, 600, 276, 3697, 839, 2138, 4504, 695, 3530, 3938, \
		1191, 1661, 617, 53, 86, 687, 1062, 1897, 2848, 1353, \
		3863, 1264, 3525, 320, 1989, 3510, 3440, 2552, 1044, 3542, \
		2774, 226, 2792, 264, 2617, 122, 1494, 3336, 876, 3812, \
		2863, 1923, 2003, 127, 539, 2305, 1094, 1263, 3891, 1354, \
		2698, 1702, 3945, 2236, 982, 837, 3699, 3613, 3577, 1696, \
		304, 4306, 2345, 1001, 253, 1935, 3034, 1432, 1872, 3739, \
		1463, 324, 1569, 3396, 274, 602, 4257, 2717, 854, 248, \
		296, 2938, 2416, 1059, 655, 4585, 3674, 1670, 2971, 3194, \
		1425, 4460, 672, 3712, 4158, 1712, 3461, 1164, 2862, 3872, \
		877, 435, 1591, 944, 3941, 2392, 2167, 4556, 4457, 3660, \
		2958, 859, 1768, 110, 1035, 3719, 2705, 1947, 1375, 2827, \
		3954, 3502, 1288, 1618, 1399, 365, 1293, 1445, 2233, 1606, \
		3055, 985, 2222, 96, 4450, 952, 1716, 2731, 2520, 398, \
		1842, 3618, 3312, 1032, 4552, 1330, 239, 594, 3288, 2769, \
		2207, 3359, 1683, 2899, 2868, 1534, 3608, 735, 893, 684, \
		2454, 2462, 3386, 2674, 4096, 1869, 2925, 4217, 1478, 3012, \
		1831, 1462, 3842, 1873, 758, 4344, 1550, 75, 1632, 791, \
		2824, 3066, 1209, 1181, 3151, 2330, 4204, 1345, 4152, 1437, \
		2686, 2704, 3796, 1036, 457, 1254, 718, 780, 4157, 3818, \
		673, 4122, 2993, 1539, 1584, 4508, 57, 1030, 3314, 1415, \
		4070, 3612, 3855, 838, 3908, 277, 579, 563, 2500, 526, \
		454, 1121, 3344, 3564, 3621, 4104, 2032, 4181, 450, 2495, \
		4031, 2933, 1777, 1589, 437, 3156, 1669, 3825, 4586, 1149, \
		2041, 1640, 4108, 2487, 3483, 3916, 471, 3026, 3092, 4166, \
		2957, 3802, 4458, 1427, 988, 3334, 1496, 3077, 1855, 3046, \
		337, 2215, 463, 2999, 783, 651, 3402, 2710, 2465, 3488, \
		1155, 140, 815, 2876, 2334, 524, 2502, 355, 3599, 2183, \
		3295, 3454, 587, 2398, 1243, 2374, 4130, 4083, 3116, 4103, \
		3687, 3565, 3311, 3770, 1843, 2538, 1811, 3576, 3854, 3700, \
		4071, 2354, 734, 3755, 1535, 811, 2370, 3163, 934, 1258, \
		4169, 2182, 3633, 356, 3949, 2164, 1074, 3251, 4568, 3966, \
		2245, 743, 4339, 1007, 2547, 201, 251, 1003, 2045, 4500, \
		4541, 3374, 186, 1695, 3853, 3614, 1812, 158, 4361, 2086, \
		2291, 1313, 2583, 184, 3376, 3310, 3620, 3688, 3345, 4533, \
		887, 1278, 823, 40, 4480, 794, 873, 4135, 3472, 1467, \
		4101, 3118, 4146, 2038, 708, 1050, 1128, 3429, 2773, 3882, \
		1045, 3464, 2745, 2350, 2069, 4195, 1326, 430, 2174, 763, \
		3937, 3903, 696, 2695, 2179, 319, 3889, 1265, 331, 1370, \
		1240, 3275, 2197, 2426, 996, 1341, 4527, 584, 362, 2404, \
		3439, 3886, 1990, 4396, 72, 567, 3420, 2116, 1287, 3790, \
		3955, 1982, 2594, 2285, 726, 3329, 2484, 1176, 2630, 4332, \
		3293, 2185, 1154, 3642, 2466, 349, 1457, 3915, 3667, 2488, \
		1733, 1624, 2151, 33, 136, 3084, 1699, 1381, 1466, 3553, \
		4136, 901, 1248, 3107, 4301, 4036, 2744, 3540, 1046, 1163, \
		3815, 1713, 2259, 389, 2838, 360, 586, 3630, 3296, 3148, \
		3258, 680, 3992, 3141, 1224, 417, 1967, 1529, 4, 919, \
		2551, 3885, 3511, 2405, 1102, 950, 4452, 1365, 448, 4183, \
		4310, 2772, 3544, 1129, 777, 1730, 2327, 1908, 2834, 4375, \
		2115, 3505, 568, 2591, 481, 2106, 1097, 1962, 3301, 2975, \
		3381, 863, 1441, 1335, 190, 1972, 4471, 1930, 2709, 3645, \
		652, 690, 2931, 4033, 273, 3838, 1570, 2447, 4020, 1926, \
		2636, 1789, 2158, 2389, 2673, 3749, 2463, 2712, 1524, 862, \
		3411, 2976, 3019, 2146, 3309, 3567, 185, 3580, 4542, 1627, \
		2367, 4497, 2025, 1449, 1144, 4175, 4415, 1680, 1230, 2910, \
		1228, 1682, 3760, 2208, 346, 4420, 4565, 2572, 2251, 3220, \
		1070, 1317, 2192, 964, 2218, 4532, 3563, 3689, 1122, 1485, \
		1796, 1912, 4580, 4133, 875, 3874, 1495, 3656, 989, 423, \
		1636, 2483, 3496, 727, 701, 1067, 4261, 308, 105, 1950, \
		3286, 596, 4046, 2575, 4191, 2751, 1414, 3703, 1031, 3769, \
		3619, 3566, 3377, 2147, 4094, 2676, 2474, 1088, 802, 2974, \
		3413, 1963, 4225, 806, 3147, 3453, 3631, 2184, 3491, 4333, \
		4199, 3183, 2768, 3763, 595, 3321, 1951, 1391, 4061, 4349, \
		382, 1835, 2299, 1021, 2009, 2196, 3520, 1241, 2400, 3961, \
		1979, 212, 2997, 465, 267, 395, 1056, 4162, 3101, 38, \
		825, 4353, 679, 3451, 3149, 1183, 4402, 1686, 2570, 4567, \
		3594, 1075, 65, 1740, 386, 867, 440, 3197, 235, 1460, \
		1833, 384, 1742, 699, 729, 1894, 2414, 2940, 666, 2885, \
		2866, 2901, 300, 3082, 138, 1157, 4144, 3120, 2715, 4259, \
		1069, 3352, 2252, 4484, 2643, 798, 1765, 1527, 1969, 1891, \
		2891, 1117, 3126, 4269, 2604, 179, 2322, 1420, 43, 314, \
		2639, 1307, 531, 234, 3244, 441, 1424, 3822, 2972, 804, \
		4227, 2189, 632, 661, 4058, 4359, 160, 2767, 3290, 4200, \
		1755, 218, 1720, 1839, 1276, 889, 522, 2336, 1771, 2268, \
		1387, 4549, 113, 4149, 1396, 770, 132, 933, 3604, 2371, \
		61, 1556, 747, 2720, 1668, 3676, 438, 869, 1906, 2329, \
		3727, 1182, 3257, 3452, 3297, 807, 2359, 606, 2565, 1223, \
		3448, 3993, 547, 4215, 2927, 2458, 4241, 1106, 677, 4355, \
		1349, 2759, 851, 750, 4268, 3209, 1118, 1039, 491, 1522, \
		2714, 3224, 4145, 3550, 4102, 3623, 4084, 373, 848, 1580, \
		152, 31, 2153, 4300, 3468, 1249, 2796, 2578, 2422, 37, \
		3263, 4163, 1517, 4468, 716, 1256, 936, 1515, 4165, 3663, \
		3027, 2765, 162, 4423, 3080, 302, 1698, 3476, 137, 3228, \
		301, 3087, 4424, 1854, 3654, 1497, 427, 1492, 124, 405, \
		148, 4099, 1469, 1878, 1208, 3730, 2825, 1377, 2623, 1147, \
		4588, 2297, 1837, 1722, 835, 984, 3781, 1607, 1597, 883, \
		3008, 2726, 2230, 3987, 336, 3652, 1856, 1561, 2079, 246, \
		856, 4517, 4026, 2093, 582, 4529, 1431, 3845, 1936, 2082, \
		2734, 1547, 2513, 2764, 3091, 3664, 472, 2379, 4266, 752, \
		1196, 2145, 3379, 2977, 2278, 2227, 4524, 2853, 1830, 3742, \
		1479, 4440, 2725, 3051, 884, 84, 55, 4510, 2813, 913, \
		4155, 782, 3648, 464, 3269, 213, 2983, 1538, 3709, 4123, \
		1977, 3963, 1751, 4233, 644, 2357, 809, 1537, 2995, 214, \
		2478, 664, 2942, 2277, 3018, 3380, 3412, 3302, 803, 3193, \
		3823, 1671, 2788, 2857, 1111, 2441, 48, 1219, 2225, 2280, \
		476, 4515, 858, 3801, 3661, 4167, 1260, 2109, 199, 2549, \
		921, 484, 1809, 2540, 1656, 2816, 516, 223, 2276, 2979, \
		665, 3234, 2415, 3830, 297, 2050, 1644, 1776, 3680, 4032, \
		3399, 691, 3975, 2457, 3137, 4216, 3745, 1870, 1434, 916, \
		767, 1621, 1804, 118, 2663, 1957, 2742, 4038, 342, 175, \
		1227, 3362, 1231, 832, 1693, 188, 1337, 2022, 2048, 299, \
		3230, 2867, 3758, 1684, 4404, 1013, 739, 4288, 1506, 1116, \
		3211, 1892, 731, 647, 1921, 2865, 3232, 667, 379, 4000, \
		2172, 432, 775, 1131, 2333, 3638, 816, 1861, 2691, 4372, \
		1412, 2753, 1533, 3757, 2900, 3231, 2886, 1922, 3871, 3813, \
		1165, 2265, 2102, 1110, 2968, 2789, 2076, 1829, 3014, 4525, \
		1343, 4206, 1352, 3893, 1898, 221, 518, 2239, 2524, 21, \
		627, 1602, 359, 3457, 390, 1410, 4374, 3423, 1909, 4330, \
		2632, 1851, 2061, 3953, 3792, 1376, 3065, 3731, 792, 4482, \
		2254, 1418, 2324, 2491, 515, 2946, 1657, 912, 3003, 4511, \
		2035, 116, 1806, 2562, 4223, 1965, 419, 4275, 1954, 1201, \
		2434, 182, 2585, 4189, 2577, 3105, 1250, 2680, 263, 3879, \
		227, 2075, 2856, 2969, 1672, 4474, 1782, 2030, 4106, 1642, \
		2052, 976, 1016, 4434, 4386, 2274, 225, 3881, 3543, 3430, \
		4311, 2206, 3762, 3289, 3184, 161, 3090, 3028, 2514, 1737, \
		1578, 850, 3130, 1350, 4208, 1311, 2293, 1532, 2870, 1413, \
		3316, 4192, 820, 974, 2054, 2349, 3539, 3465, 4037, 2915, \
		1958, 1761, 2531, 166, 1513, 938, 1546, 3031, 2083, 2519, \
		3774, 1717, 1901, 4522, 2229, 3050, 3009, 4441, 494, 845, \
		1667, 3158, 748, 853, 3834, 4258, 3223, 3121, 1523, 3384, \
		2464, 3644, 3403, 1931, 2527, 1946, 3795, 3720, 2687, 4399, \
		2621, 1379, 1701, 3861, 1355, 2178, 3528, 697, 1744, 4371, \
		2873, 1862, 70, 4398, 2703, 3721, 1438, 2263, 1167, 1758, \
		262, 2794, 1251, 12, 2473, 3306, 4095, 3748, 3387, 2390, \
		3943, 1704, 1727, 721, 2589, 570, 1199, 1956, 2917, 119, \
		1187, 1883, 2065, 4572, 4337, 745, 1558, 2510, 4347, 4063, \
		4139, 2098, 4078, 2162, 3951, 2063, 1885, 797, 3217, 4485, \
		1269, 1306, 3201, 315, 1788, 3391, 1927, 1675, 1850, 2831, \
		4331, 3493, 1177, 1291, 367, 1785, 4173, 1146, 3063, 1378, \
		2701, 4400, 1185, 121, 3877, 265, 467, 1600, 629, 1320, \
		4380, 635, 1091, 1998, 4112, 414, 178, 3207, 4270, 704, \
		1082, 196, 537, 129, 3921, 550, 2284, 3499, 1983, 480, \
		3418, 569, 2667, 722, 4465, 4188, 2799, 183, 3569, 1314, \
		2396, 589, 2421, 3104, 2797, 4190, 3318, 4047, 2250, 3354, \
		4566, 3253, 1687, 2544, 93, 1222, 3143, 607, 4222, 2808, \
		1807, 486, 906, 411, 79, 1780, 4476, 444, 1043, 3884, \
		3441, 920, 2952, 200, 3587, 1008, 92, 2568, 1688, 4293, \
		1655, 2948, 1810, 3616, 1844, 2243, 3968, 4437, 4562, 165, \
		2739, 1762, 639, 1945, 2707, 1932, 20, 2843, 2240, 1054, \
		397, 3773, 2732, 2084, 4363, 1801, 1736, 2763, 3029, 1548, \
		4346, 2654, 1559, 1858, 2073, 229, 1019, 2301, 354, 3635, \
		525, 3693, 564, 1553, 4010, 4030, 3682, 451, 1941, 514, \
		2818, 2325, 1732, 3482, 3668, 4109, 1175, 3495, 3330, 1637, \
		4393, 4056, 663, 2981, 215, 1170, 1087, 3305, 2677, 13, \
		4264, 2381, 4247, 4418, 348, 3487, 3643, 2711, 3385, 3750, \
		2455, 3977, 4240, 3136, 2928, 3976, 2461, 3751, 685, 88, \
		620, 208, 3933, 4019, 3394, 1571, 1213, 1283, 498, 47, \
		2966, 1112, 4559, 1482, 292, 2320, 181, 2801, 1202, 2203, \
		4015, 512, 1943, 641, 995, 3518, 2198, 786, 36, 3103, \
		2579, 590, 1710, 4160, 1058, 3829, 2939, 3235, 1895, 1064, \
		4273, 421, 991, 959, 1297, 1101, 3438, 3512, 363, 1401, \
		3960, 3273, 1242, 3628, 588, 2581, 1315, 1072, 2166, 3806, \
		3942, 2672, 3388, 2159, 577, 279, 1821, 830, 1233, 4246, \
		2470, 4265, 3024, 473, 4253, 941, 4129, 3626, 1244, 60, \
		3162, 3605, 812, 4496, 3371, 1628, 2312, 1575, 68, 1864, \
		1542, 605, 3145, 808, 2986, 645, 733, 3610, 4072, 1747, \
		2068, 3538, 2746, 2055, 624, 1000, 3849, 4307, 4446, 543, \
		2018, 1238, 1372, 108, 1770, 3174, 523, 3637, 2877, 1132, \
		4203, 3726, 3152, 1907, 3425, 1731, 2490, 2819, 1419, 3205, \
		180, 2436, 293, 204, 1815, 4367, 1648, 1613, 1574, 2365, \
		1629, 4116, 340, 4040, 1996, 1093, 3866, 540, 612, 353, \
		2504, 1020, 3279, 1836, 3060, 4589, 2, 1531, 2755, 1312, \
		3571, 2087, 4238, 3979, 4051, 725, 3498, 2595, 551, 4251, \
		475, 2962, 2226, 3017, 2978, 2943, 224, 2776, 4387, 4298, \
		2155, 4411, 1386, 3172, 1772, 2101, 2860, 1166, 2684, 1439, \
		865, 388, 3459, 1714, 954, 4068, 1417, 2821, 4483, 3219, \
		3353, 2573, 4048, 145, 284, 742, 3591, 3967, 2536, 1845, \
		1053, 2523, 2844, 519, 981, 3858, 3946, 1605, 3783, 1446, \
		3986, 3049, 2727, 4523, 3016, 2279, 2963, 1220, 95, 3779, \
		986, 1429, 4531, 3347, 965, 462, 3650, 338, 4118, 1567, \
		326, 2121, 345, 3358, 3761, 2770, 4312, 4014, 2432, 1203, \
		1919, 649, 785, 2425, 3519, 3276, 2010, 4003, 963, 3349, \
		1318, 631, 3190, 4228, 1079, 1153, 3490, 3294, 3632, 3600, \
		4170, 318, 3527, 2696, 1356, 1652, 762, 3533, 431, 2881, \
		4001, 2012, 1488, 4555, 3805, 2393, 1073, 3596, 3950, 2648, \
		4079, 576, 2388, 3389, 1790, 4410, 2271, 4299, 3109, 32, \
		3479, 1625, 4544, 4093, 3308, 3378, 3020, 1197, 572, 658, \
		4383, 3971, 4503, 3906, 840, 2130, 615, 1663, 1455, 351, \
		614, 2136, 841, 1206, 1880, 948, 1104, 4243, 173, 344, \
		2210, 327, 4578, 1914, 1286, 3504, 3421, 4376, 1610, 3928, \
		535, 198, 2954, 1261, 1096, 3416, 482, 923, 1109, 2859, \
		2266, 1773, 4077, 2650, 4140, 909, 561, 581, 3038, 4027, \
		898, 4066, 956, 4237, 2290, 3572, 4362, 2518, 2733, 3032, \
		1937, 245, 3043, 1562, 1828, 2855, 2790, 228, 2507, 1859, \
		818, 4194, 3537, 2351, 1748, 4571, 2659, 1884, 2646, 3952, \
		2829, 1852, 4426, 29, 154, 623, 2348, 2747, 975, 2781, \
		1643, 2936, 298, 2903, 2023, 4499, 3583, 1004, 4391, 1639, \
		3671, 1150, 707, 3548, 4147, 115, 2811, 4512, 4180, 3685, \
		4105, 2784, 1783, 369, 3984, 1448, 3369, 4498, 2047, 2904, \
		1338, 25, 1237, 2341, 544, 1273, 401, 1794, 1487, 2170, \
		4002, 2195, 3277, 1022, 4408, 1792, 403, 126, 3869, 1924, \
		4022, 1173, 4111, 2608, 1092, 2307, 4041, 670, 4462, 4054, \
		4395, 3509, 3887, 321, 1384, 4413, 4177, 479, 2593, 3500, \
		3956, 211, 3271, 3962, 2991, 4124, 4326, 714, 4470, 3406, \
		191, 1890, 3213, 1528, 3445, 418, 2806, 4224, 3300, 3414, \
		1098, 260, 1760, 2741, 2916, 2664, 1200, 2803, 4276, 1390, \
		3285, 3322, 106, 1374, 3794, 2706, 2528, 640, 2429, 513, \
		2493, 452, 528, 244, 2081, 3033, 3846, 254, 19, 2526, \
		2708, 3404, 4472, 1674, 2635, 3392, 4021, 2002, 3870, 2864, \
		2887, 648, 2201, 1204, 843, 496, 1285, 2118, 4579, 3340, \
		1797, 4329, 2833, 3424, 2328, 3153, 870, 1888, 193, 4521, \
		2729, 1718, 220, 2847, 3894, 1063, 2413, 3236, 730, 2890, \
		3212, 1970, 192, 1904, 871, 796, 2645, 2064, 2660, 1188, \
		947, 2127, 1207, 3068, 1470, 1368, 333, 757, 3738, 3843, \
		1433, 2924, 3746, 4097, 150, 1582, 1541, 2362, 69, 2690, \
		2874, 817, 2072, 2508, 1560, 3045, 3653, 3078, 4425, 2060, \
		2830, 2633, 1676, 554, 1126, 1052, 2242, 2537, 3617, 3771, \
		399, 1275, 3178, 1721, 3059, 2298, 3280, 383, 3241, 1461, \
		3741, 3013, 2854, 2077, 1563, 1501, 4430, 1407, 1302, 829, \
		2385, 280, 1161, 1048, 710, 4366, 2317, 205, 157, 3575, \
		3615, 2539, 2949, 485, 2561, 2809, 117, 2919, 1622, 1735, \
		2516, 4364, 712, 4328, 1911, 3341, 1486, 2014, 402, 2006, \
		4409, 2157, 3390, 2637, 316, 4172, 2626, 368, 2029, 2785, \
		4475, 2556, 80, 1588, 3679, 2934, 1645, 4076, 2100, 2267, \
		3173, 2337, 109, 3799, 860, 1526, 3215, 799, 638, 2530, \
		2740, 1959, 261, 2682, 1168, 217, 3181, 4201, 1134, 4232, \
		2989, 3964, 4570, 2067, 2352, 4073, 4370, 2693, 698, 3239, \
		385, 3248, 66, 1577, 2762, 2515, 1802, 1623, 3481, 2489, \
		2326, 3426, 778, 720, 2669, 1705, 1360, 1691, 834, 3058, \
		1838, 3179, 219, 1900, 2730, 3775, 953, 2258, 3460, 3816, \
		4159, 2419, 591, 4212, 3924, 1359, 1726, 2670, 3944, 3860, \
		2699, 1380, 3475, 3085, 303, 3852, 3578, 187, 2907, 833, \
		1724, 1361, 4292, 2543, 2569, 3254, 4403, 2898, 3759, 3360, \
		1229, 3364, 4416, 4249, 553, 1849, 2634, 1928, 4473, 2787, \
		2970, 3824, 3675, 3157, 2721, 846, 375, 1454, 2134, 616, \
		3900, 1192, 559, 911, 2815, 2947, 2541, 4294, 761, 2176, \
		1357, 3926, 1612, 2315, 4368, 4075, 1775, 2935, 2051, 2782, \
		4107, 3670, 2042, 4392, 2482, 3331, 424, 930, 790, 3733, \
		76, 4115, 2311, 2366, 3372, 4543, 2150, 3480, 1734, 1803, \
		2920, 768, 1398, 3788, 1289, 1179, 1211, 1573, 2314, 1649, \
		3927, 2113, 4377, 1596, 3054, 3782, 2234, 3947, 358, 2840, \
		628, 2614, 468, 882, 3053, 1608, 4378, 1322, 4127, 943, \
		3809, 436, 3678, 1778, 81, 4536, 4507, 3707, 1540, 1866, \
		151, 3112, 849, 2761, 1738, 67, 2364, 2313, 1614, 1212, \
		2446, 3395, 3839, 325, 2212, 4119, 927, 1500, 1827, 2078, \
		3044, 1857, 2509, 2655, 746, 3160, 62, 4009, 2498, 565, \
		74, 3735, 4345, 2512, 3030, 2735, 939, 4255, 604, 2361, \
		1865, 1583, 3708, 2994, 2984, 810, 3607, 3756, 2869, 2754, \
		2294, 3, 3444, 1968, 3214, 1766, 861, 3383, 2713, 3122, \
		492, 4443, 4186, 4467, 3099, 4164, 3094, 937, 2737, 167, \
		1504, 4290, 1363, 4454, 1115, 2893, 4289, 1511, 168, 4429, \
		1826, 1564, 928, 426, 3076, 3655, 3335, 3875, 123, 3074, \
		428, 1328, 4554, 2169, 2013, 1795, 3342, 1123, 291, 2438, \
		4560, 4439, 3011, 3743, 4218, 507, 968, 489, 1041, 446, \
		1367, 1877, 3069, 4100, 3552, 3473, 1382, 323, 3841, 3740, \
		1832, 3242, 236, 3914, 3485, 350, 2133, 1664, 376, 4044, \
		598, 1143, 3368, 2026, 3985, 2232, 3784, 1294, 1139, 1334, \
		3409, 864, 2262, 2685, 3722, 4153, 915, 2923, 1871, 3844, \
		3035, 4530, 2220, 987, 3658, 4459, 3821, 3195, 442, 4478, \
		42, 3204, 2323, 2820, 2255, 4069, 3702, 3315, 2752, 2871, \
		4373, 2836, 391, 1301, 1824, 4431, 232, 533, 3930, 3959, \
		2402, 364, 3787, 1619, 769, 3167, 4150, 1347, 4357, 4060, \
		3284, 1952, 4277, 4548, 3171, 2269, 4412, 1987, 322, 1465, \
		3474, 1700, 2700, 2622, 3064, 2826, 3793, 1948, 107, 2339, \
		1239, 3522, 332, 1876, 1471, 447, 3434, 4453, 1509, 4291, \
		1690, 1725, 1706, 3925, 1651, 2177, 2697, 3862, 3892, 2849, \
		4207, 2758, 3131, 4356, 1394, 4151, 3724, 4205, 2851, 4526, \
		3516, 997, 24, 2021, 2905, 189, 3408, 1442, 1140, 3912, \
		238, 3766, 4553, 1490, 429, 3535, 4196, 4324, 4126, 1594, \
		4379, 2612, 630, 2191, 3350, 1071, 2395, 2582, 3570, 2292, \
		2756, 4209, 242, 530, 3200, 2640, 1270, 3996, 828, 1823, \
		1408, 392, 258, 1100, 2407, 960, 1138, 1444, 3785, 366, \
		2628, 1178, 1617, 3789, 3503, 2117, 1915, 497, 2444, 1214, \
		1025, 972, 822, 3560, 888, 3177, 1840, 400, 2016, 545, \
		3995, 1305, 2641, 4486, 4283, 330, 3524, 3890, 3864, 1095, \
		2108, 2955, 4168, 3602, 935, 3096, 717, 3716, 458, 11, \
		2679, 2795, 3106, 3469, 902, 1028, 59, 2373, 3627, 2399, \
		3274, 3521, 1371, 2340, 2019, 26, 171, 4245, 2383, 831, \
		2909, 3363, 1681, 3361, 2911, 176, 416, 3447, 3142, 2566, \
		94, 2224, 2964, 49, 1011, 4406, 1024, 1282, 2445, 1572, \
		1615, 1180, 3729, 3067, 1879, 2128, 842, 1918, 2202, 2433, \
		2802, 1955, 2665, 571, 2144, 3021, 753, 4318, 558, 1660, \
		3901, 3939, 946, 1882, 2661, 120, 2619, 4401, 3256, 3150, \
		3728, 1210, 1616, 1290, 2629, 3494, 2485, 4110, 2000, 4023, \
		1086, 2476, 216, 1757, 2683, 2264, 2861, 3814, 3462, 1047, \
		1819, 281, 408, 4143, 3226, 139, 3641, 3489, 2186, 1080, \
		706, 2040, 3672, 4587, 3062, 2624, 4174, 3367, 1450, 599, \
		3911, 1333, 1443, 1295, 961, 4005, 4231, 1753, 4202, 2332, \
		2878, 776, 3428, 3545, 1051, 1847, 555, 290, 1484, 3343, \
		3690, 455, 1038, 3125, 3210, 2892, 1507, 4455, 4558, 2440, \
		2967, 2858, 2103, 924, 676, 3134, 4242, 2125, 949, 3437, \
		2406, 1298, 259, 1961, 3415, 2107, 1262, 3865, 2306, 1997, \
		2609, 636, 801, 3304, 2475, 1171, 4024, 4519, 195, 2601, \
		705, 1152, 2187, 4229, 4007, 64, 3250, 3595, 2165, 2394, \
		1316, 3351, 3221, 4260, 3326, 702, 4272, 2412, 1896, 3895, \
		688, 654, 3828, 2417, 4161, 3265, 396, 2522, 2241, 1846, \
		1127, 3546, 709, 1818, 1162, 3463, 3541, 3883, 2553, 445, \
		1473, 490, 3124, 1119, 456, 3718, 3797, 111, 4551, 3768, \
		3313, 3704, 58, 1246, 903, 971, 1281, 1215, 4407, 2008, \
		3278, 2300, 2505, 230, 4433, 2779, 977, 738, 2896, 4405, \
		1217, 50, 91, 2546, 3588, 4340, 4390, 2044, 3584, 252, \
		3848, 2346, 625, 23, 1340, 3517, 2427, 642, 4235, 958, \
		2409, 422, 3333, 3657, 1428, 2221, 3780, 3056, 836, 3857, \
		2237, 520, 891, 737, 1015, 2780, 2053, 2748, 821, 1280, \
		1026, 904, 488, 1475, 508, 461, 2217, 3348, 2193, 4004, \
		1137, 1296, 2408, 992, 4236, 2089, 4067, 2257, 1715, 3776, \
		4451, 3436, 1103, 2126, 1881, 1189, 3940, 3808, 1592, 4128, \
		2376, 4254, 1545, 2736, 1514, 3095, 1257, 3603, 3164, 133, \
		789, 1634, 425, 1499, 1565, 4120, 675, 1108, 2104, 483, \
		2951, 2550, 3442, 5, 766, 2922, 1435, 4154, 3002, 2814, \
		1658, 560, 2096, 4141, 410, 2559, 487, 970, 1027, 1247, \
		3470, 4137, 4065, 2091, 4028, 4012, 4314, 683, 3753, 736, \
		979, 521, 3176, 1277, 3561, 4534, 83, 3007, 3052, 1598, \
		469, 3918, 773, 434, 3811, 3873, 3337, 4134, 3555, 795, \
		1887, 1905, 3154, 439, 3246, 387, 2261, 1440, 3410, 3382, \
		1525, 1767, 3800, 2959, 4516, 3041, 247, 3833, 2718, 749, \
		3129, 2760, 1579, 3113, 374, 1666, 2722, 495, 1917, 1205, \
		2129, 2137, 3907, 3698, 3856, 983, 3057, 1723, 1692, 2908, \
		1232, 2384, 1822, 1303, 3997, 4352, 3261, 39, 3559, 1279, \
		973, 2749, 4193, 2071, 1860, 2875, 3639, 141, 4495, 2369, \
		3606, 1536, 2985, 2358, 3146, 3298, 4226, 3192, 2973, 3303, \
		1089, 637, 1764, 3216, 2644, 1886, 872, 3556, 4481, 2823, \
		3732, 1633, 931, 134, 35, 2424, 2199, 650, 3647, 3000, \
		4156, 3714, 719, 1729, 3427, 1130, 2879, 433, 879, 3919, \
		131, 3166, 1397, 1620, 2921, 917, 6, 3936, 3532, 2175, \
		1653, 4295, 4343, 3737, 1874, 334, 3989, 4317, 1195, 3022, \
		4267, 3128, 852, 2719, 3159, 1557, 2656, 4338, 3590, 2246, \
		285, 4287, 2895, 1014, 978, 892, 3754, 3609, 2355, 646, \
		2889, 1893, 3237, 700, 3328, 3497, 2286, 4052, 4464, 2588, \
		2668, 1728, 779, 3715, 1255, 3097, 4469, 1974, 4327, 1799, \
		4365, 1817, 1049, 3547, 2039, 1151, 1081, 2602, 4271, 1066, \
		3327, 728, 3238, 1743, 2694, 3529, 3904, 4505, 4538, 3974, \
		2930, 3400, 653, 1061, 3896, 87, 2453, 3752, 894, 4315, \
		3991, 3450, 3259, 4354, 3133, 1107, 925, 4121, 3711, 3819, \
		4461, 1994, 4042, 378, 2884, 3233, 2941, 2980, 2479, 4057, \
		3188, 633, 4382, 2142, 573, 4584, 3827, 1060, 689, 3401, \
		3646, 784, 2200, 1920, 2888, 732, 2356, 2987, 4234, 994, \
		2428, 1944, 2529, 1763, 800, 1090, 2610, 4381, 660, 3189, \
		2190, 1319, 2613, 1601, 2841, 22, 999, 2347, 2056, 155, \
		207, 2451, 89, 52, 3899, 1662, 2135, 2131, 352, 2303, \
		541, 4448, 98, 4221, 2564, 3144, 2360, 1543, 4256, 3836, \
		275, 3910, 1142, 1451, 4045, 3320, 3287, 3764, 240, 4211, \
		1709, 2420, 2580, 2397, 3629, 3455, 361, 3514, 4528, 3037, \
		2094, 562, 3695, 278, 2387, 2160, 4080, 4583, 657, 2143, \
		1198, 2666, 2590, 3419, 3506, 73, 1552, 2499, 3694, 580, \
		2095, 910, 1659, 1193, 4319, 289, 1125, 1848, 1677, 4250, \
		2283, 2596, 3922, 4214, 3139, 3994, 1272, 2017, 2342, 4447, \
		611, 2304, 3867, 128, 2599, 197, 2111, 3929, 1404, 233, \
		3199, 1308, 243, 1939, 453, 3692, 2501, 3636, 2335, 3175, \
		890, 980, 2238, 2845, 222, 2945, 2817, 2492, 1942, 2430, \
		4016, 9, 460, 967, 1476, 4219, 100, 501, 4089, 4088, \
		504, 101, 46, 2443, 1284, 1916, 844, 2723, 4442, 1521, \
		3123, 1040, 1474, 969, 905, 2560, 1808, 2950, 922, 2105, \
		3417, 2592, 1984, 4178, 4514, 2961, 2281, 4252, 2378, 3025, \
		3665, 3917, 881, 1599, 2615, 266, 3268, 2998, 3649, 2216, \
		966, 509, 10, 1253, 3717, 1037, 1120, 3691, 527, 1940, \
		2494, 3683, 4182, 3433, 1366, 1472, 1042, 2554, 4477, 1423, \
		3196, 3245, 868, 3155, 3677, 1590, 3810, 878, 774, 2880, \
		2173, 3534, 1327, 1491, 3075, 1498, 929, 1635, 3332, 990, \
		2410, 4274, 2805, 1966, 3446, 1225, 177, 2606, 4113, 78, \
		2558, 907, 4142, 1159, 282, 147, 3072, 125, 2005, 1793, \
		2015, 1274, 1841, 3772, 2521, 1055, 3266, 268, 257, 1300, \
		1409, 2837, 3458, 2260, 866, 3247, 1741, 3240, 1834, 3281, \
		4350, 3999, 2883, 668, 4043, 1453, 1665, 847, 3114, 4085, \
		4492, 3983, 2028, 1784, 2627, 1292, 3786, 1400, 2403, 3513, \
		585, 3456, 2839, 1603, 3948, 3598, 3634, 2503, 2302, 613, \
		2132, 1456, 3486, 2467, 4419, 3357, 2209, 2122, 174, 2913, \
		4039, 2309, 4117, 2214, 3651, 3047, 3988, 756, 1875, 1369, \
		3523, 1266, 4284, 4577, 2120, 2211, 1568, 3840, 1464, 1383, \
		1988, 3888, 3526, 2180, 4171, 1787, 2638, 3202, 44, 103, \
		310, 311, 104, 3324, 4262, 15, 4305, 3851, 1697, 3086, \
		3081, 3229, 2902, 2049, 2937, 3831, 249, 203, 2319, 2437, \
		1483, 1124, 556, 4320, 4575, 4286, 741, 2247, 146, 407, \
		1160, 1820, 2386, 578, 3696, 3909, 601, 3837, 3397, 4034, \
		4303, 17, 256, 394, 3267, 466, 2616, 3878, 2793, 2681, \
		1759, 1960, 1099, 1299, 393, 269, 18, 1934, 3847, 1002, \
		3585, 202, 295, 3832, 855, 3042, 2080, 1938, 529, 1309, \
		4210, 593, 3765, 1331, 3913, 1459, 3243, 3198, 532, 1405, \
		4432, 1018, 2506, 2074, 2791, 3880, 2775, 2275, 2944, 517, \
		2846, 1899, 1719, 3180, 1756, 1169, 2477, 2982, 2996, 3270, \
		1980, 3957, 3932, 2450, 621, 156, 1814, 2318, 294, 250, \
		3586, 2548, 2953, 2110, 536, 2600, 1083, 4520, 1903, 1889, \
		1971, 3407, 1336, 2906, 1694, 3579, 3375, 3568, 2584, 2800, \
		2435, 2321, 3206, 2605, 415, 1226, 2912, 343, 2123, 4244, \
		1235, 27, 4428, 1503, 1512, 2738, 2532, 4563, 4422, 3089, \
		2766, 3185, 4360, 3574, 1813, 206, 622, 2057, 30, 3111, \
		1581, 1867, 4098, 3071, 406, 283, 2248, 4049, 3981, 4494, \
		814, 3640, 1156, 3227, 3083, 3477, 34, 788, 932, 3165, \
		771, 3920, 2598, 538, 3868, 2004, 404, 3073, 1493, 3876, \
		2618, 1186, 2662, 2918, 1805, 2810, 2036, 4148, 3169, 4550, \
		1034, 3798, 1769, 2338, 1373, 1949, 3323, 309, 312, 45, \
		500, 505, 4220, 609, 4449, 3778, 2223, 1221, 2567, 2545, \
		1009, 51, 619, 2452, 686, 3897, 54, 3006, 885, 4535, \
		1587, 1779, 2557, 412, 4114, 1631, 3734, 1551, 566, 3507, \
		4397, 2689, 1863, 2363, 1576, 1739, 3249, 1076, 4008, 1555, \
		3161, 2372, 1245, 1029, 3705, 4509, 3005, 85, 3898, 618, \
		90, 1010, 1218, 2965, 2442, 499, 102, 313, 3203, 1421, \
		4479, 3558, 824, 3262, 3102, 2423, 787, 135, 3478, 2152, \
		3110, 153, 2058, 4427, 170, 1236, 2020, 1339, 998, 626, \
		2842, 2525, 1933, 255, 270, 4304, 306, 4263, 2472, 2678, \
		1252, 459, 510, 4017, 3935, 765, 918, 3443, 1530, 2295, \
		4590};

// ========================================
//   Helper functions
// ========================================

/*
 *   poly_init:
 *   	Initialize a polynomial with a constant.
 */

Poly* poly_init(int value) {
	Poly* p = malloc(sizeof(Poly));
	p->degree = 0;
	p->coefficients = malloc(sizeof(int));
	p->coefficients[0] = value;
	p->allocated_num_coeffs = 1;
	return p;
}

/*
 *   poly_gen:
 *   	Generates a random polynomial of degree n.
 */
Poly* poly_gen() {
	Poly* p = poly_init(1);
	poly_assert_degree(p, INPUT_LENGTH);
	p->degree = INPUT_LENGTH;
	for (int i = 0; i <= p->degree; i++) {
		p->coefficients[i] = rand() % INPUT_BASE;
	}
	return p;
}

void poly_free(Poly* p) {
	free(p->coefficients);
	free(p);
}

/*
 *   poly_assign:
 *   	Copies a polynomial to destination.
 */

void poly_assign(Poly* dst, const Poly* src) {
	poly_assert_degree(dst, src->degree);

	for (int i = 0; i <= src->degree; i++) {
		dst->coefficients[i] = src->coefficients[i];
	}

	dst->degree = src->degree;
}

/*
 *   poly_read:
 *   	Reads input string as coefficients of a polynomial.
 *   	In this case, each coefficient is a digit between 0-9.
 *   	Note that throughout the computation, the polynomials are in ascending order.
 */
void poly_read(Poly* dst, const char* src) {
	int n = strlen(src);
	poly_assert_degree(dst, n-1);

	for (int i = 0; i < n; i++) {
		dst->coefficients[i] = char2int(src[i]);
	}
	dst->degree = n-1;
}

/*
 *   poly_assert_degree:
 *   	Makes sure p has enough space to hold target polynomial.
 *   	If not, reallocate more space with values initiated to 0.
 */
void poly_assert_degree(Poly* p, const int degree_needed) {
	int alloc_before = p->allocated_num_coeffs;
	if (alloc_before < degree_needed + 1) {
		p->coefficients = realloc(p->coefficients, (degree_needed + 1) * sizeof(int));
		for (int i = p->degree + 1; i <= degree_needed; i++) {
			p->coefficients[i] = 0;
		}
		p->allocated_num_coeffs = degree_needed + 1;
	}
}

void poly_assert_equal_degree(Poly* r, Poly* p, Poly* q) {
	int res_degree = MAX(p->degree, q->degree);
	poly_assert_degree(r, res_degree);
	poly_assert_degree(p, res_degree);
	poly_assert_degree(q, res_degree);
}

/*
 *   poly_slice:
 *   	Slice the given polynomial.
 */
void poly_slice(Poly* dst, Poly* src, int begin, int degree) {
	Poly* result = poly_init(0);
	if (degree > src->degree) {
		poly_assign(dst, src);
	} else {
		poly_assert_degree(result, degree);
		int i;
		result->degree = degree;
		for (i = 0; i <= result->degree; i++) {
			result->coefficients[i] = src->coefficients[i+begin];
		}
		poly_assign(dst, result);
	}
	poly_free(result);
}

void poly_print(const Poly* p) {
	// constant term
	printf("%d ", p->coefficients[0]);
	// general terms
	for (int i = 1; i <= p->degree; i++) {
		if (p->coefficients[i] >= 0) {
			printf("+ %d", p->coefficients[i]);
		} else {
			printf("- %d", p->coefficients[i]*(-1));
		}
		printf("*x^%d ", i);
	}
	printf("\n");
}

/*
 *   poly_rshift/poly_lshift:
 *   	Shifts the polynomial to perform multiplication/division of x^k.
 */
void poly_rshift(Poly* dst, Poly* src, int bias) {
	Poly* result = poly_init(0);
	int res_degree = src->degree + bias;
	poly_assert_degree(result, res_degree);
	int i;
	for (i = 0; i < bias; i++) {
		result->coefficients[i] = 0;
	}
	for (i = 0; i + bias <= res_degree; i++) {
		result->coefficients[i + bias] = src->coefficients[i];
	}
	result->degree = res_degree;

	poly_assign(dst, result);
	poly_free(result);
}

void poly_lshift(Poly* dst, Poly* src, int bias) {
	Poly* result = poly_init(0);
	int res_degree = src->degree - bias;
	poly_assert_degree(result, res_degree);
	int i;
	for (i = 0; i <= res_degree; i++) {
		result->coefficients[i] = src->coefficients[i + bias];
	}
	result->degree = res_degree;

	poly_assign(dst, result);
	poly_free(result);
}

/*
 *   poly_compare:
 *   	Compares two polynomials.
 */
int poly_compare(const Poly* p, const Poly* q) {
	if (p->degree > q->degree) {
		return 1;
	} else if (p->degree < q->degree) {
		return -1;
	} else {
		return 0;
	}
}

// ========================================
//   Arithmetic operations
//   Results are stored in the first argument
// ========================================

void poly_add(Poly* r, Poly* p, Poly* q) {
	Poly* result = poly_init(0);
	int res_degree = MAX(p->degree, q->degree);
	poly_assert_equal_degree(result, p, q);

	for (int i = 0; i <= res_degree; i++) {
		result->coefficients[i] = p->coefficients[i] + q->coefficients[i];
	}
	result->degree = res_degree;

	poly_assign(r, result);
	poly_free(result);
}

void poly_add_mod(Poly* r, Poly* p, Poly* q) {
	Poly* result = poly_init(0);
	int res_degree = MAX(p->degree, q->degree);
	poly_assert_equal_degree(result, p, q);

	for (int i = 0; i <= res_degree; i++) {
		result->coefficients[i] = (p->coefficients[i] + q->coefficients[i]) % MODULUS;
	}
	result->degree = res_degree;

	poly_assign(r, result);
	poly_free(result);
}

void poly_add_int(Poly* r, Poly* p, const int num) {
	Poly* addend = poly_init(num);
	poly_add(r, p, addend);
	poly_free(addend);
}

void poly_add_mod_int(Poly* r, Poly* p, const int num) {
	Poly* addend = poly_init(num);
	poly_add_mod(r, p, addend);
	poly_free(addend);
}

void poly_sub(Poly* r, Poly* p, Poly* q) {
	// If deg(p)>= deg(q) => works fine
	// If deg(p) < deg(q) => need to increase p's degree
	Poly* result = poly_init(0);
	int res_degree = MAX(p->degree, q->degree);
	poly_assert_equal_degree(result, p, q);

	for (int i = 0; i <= res_degree; i++) {
		result->coefficients[i] = p->coefficients[i] - q->coefficients[i];
	}
	result->degree = res_degree;

	poly_assign(r, result);
	poly_free(result);
}

void poly_sub_mod(Poly* r, Poly* p, Poly* q) {
	Poly* result = poly_init(0);
	int res_degree = MAX(p->degree, q->degree);
	poly_assert_equal_degree(result, p, q);

	for (int i = 0; i <= res_degree; i++) {
		result->coefficients[i] = (p->coefficients[i] - q->coefficients[i]) % MODULUS;
	}
	result->degree = res_degree;

	poly_assign(r, result);
	poly_free(result);
}

void poly_sub_int(Poly* r, Poly* p, const int num) {
	Poly* subtrahend = poly_init(num);
	poly_sub(r, p, subtrahend);
	poly_free(subtrahend);
}

void poly_sub_mod_int(Poly* r, Poly* p, const int num) {
	Poly* subtrahend = poly_init(num);
	poly_sub_mod(r, p, subtrahend);
	poly_free(subtrahend);
}

void poly_mul(Poly* r, Poly* p, Poly* q) {
	// schoolbook method

	if (DEBUG) {
		printf("Entering poly_mul(): \n");
		printf("\tBefore multiplication,\n\t p is: \n\t\t"); poly_print(p);
		printf("\t q is: \n\t\t"); poly_print(q);
	}
	Poly* result = poly_init(0);
	int res_degree = p->degree + q->degree;
	poly_assert_degree(result, res_degree);

	int i, j;
	for (i = 0; i < result->allocated_num_coeffs; i++) {
		 result->coefficients[i]= 0;
	}

	for (i = 0; i <= p->degree; i++) {
		for (j = 0; j <= q->degree; j++) {
			result->coefficients[i+j] += p->coefficients[i] * q->coefficients[j];
		}
	}
	result->degree = res_degree;

	if (DEBUG) { printf("\tResult of poly_mul(): \n\t\t"); poly_print(result); }

	poly_assign(r, result);
	poly_free(result);
}

void poly_mul_mod(Poly* r, Poly* p, Poly* q) {
	// schoolbook method
	if (DEBUG) {
		printf("Entering poly_mul_mod(): \n");
		printf("\tBefore multiplication,\n\t p is: \n\t\t"); poly_print(p);
		printf("\t q is: \n\t\t"); poly_print(q);
	}
	Poly* result = poly_init(0);
	int res_degree = p->degree + q->degree;
	poly_assert_degree(result, res_degree);

	int i, j;
	for (i = 0; i < result->allocated_num_coeffs; i++) {
		 result->coefficients[i]= 0;
	}

	for (i = 0; i <= p->degree; i++) {
		for (j = 0; j <= q->degree; j++) {
			result->coefficients[i+j] += p->coefficients[i] * q->coefficients[j];
			result->coefficients[i+j] %= MODULUS;
		}
	}
	result->degree = res_degree;

	if (DEBUG) { printf("\tResult of poly_mul_mod(): \n\t\t"); poly_print(result); printf("\n"); }

	poly_assign(r, result);
	poly_free(result);
}

void poly_mul_int(Poly* r, Poly* p, const int num) {
	Poly* multiplier = poly_init(num);
	poly_mul(r, p, multiplier);
	poly_free(multiplier);
}

void poly_mul_mod_int(Poly* r, Poly* p, const int num) {
	Poly* multiplier = poly_init(num % MODULUS);
	poly_mul_mod(r, p, multiplier);
	poly_free(multiplier);
}

/*
 *  poly_div_mod_int:
 * 		Divides polynomial p by an integer under modulo n.
 * 		The division is implemented as multiplying its multiplicative inverse.
 */
void poly_div_mod_int(Poly* r, Poly* p, const int num) {
	Poly* inv_divisor = poly_init(mod_inv[num % MODULUS]);
	poly_mul_mod(r, p, inv_divisor);
	poly_free(inv_divisor);
}

// ========================================
//   Miscellaneous
// ========================================

/*
 *   poly_eval_print:
 *   	Prints coefficients of p as a large integer.
 */
void poly_eval_print(Poly* p) {
	int i, carry = 0;
	poly_assert_degree(p, p->degree+1);
	for (i = 0; i <= p->degree; i++) {
		carry = p->coefficients[i] / INPUT_BASE;
		p->coefficients[i] %= INPUT_BASE;
		p->coefficients[i+1] += carry;
	}
	for (i = p->degree; i >= 0; i--) {
		printf("%d", p->coefficients[i]);
	}
}

/*
 *   poly_eval_mod:
 *   	Evaluates p at x0 and returns.
 */
int poly_eval_mod(Poly* p , int x0) {
	int result = p->coefficients[0];
	for (int i = 1; i <= p->degree; i++) {
		result = (result * x0 + p->coefficients[i]) % MODULUS;
	}
	return result;
}

