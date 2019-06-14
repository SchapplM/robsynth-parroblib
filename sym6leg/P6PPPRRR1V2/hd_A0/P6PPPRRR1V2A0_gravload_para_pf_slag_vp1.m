% Calculate Gravitation load for parallel robot
% P6PPPRRR1V2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d4,theta1,theta2,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-16 19:43
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6PPPRRR1V2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(10,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPPRRR1V2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 19:29:02
% EndTime: 2019-05-16 19:29:09
% DurationCPUTime: 8.02s
% Computational Cost: add. (3389->426), mult. (8043->819), div. (378->3), fcn. (9779->50), ass. (0->312)
t3074 = xP(5);
t3038 = sin(t3074);
t3041 = cos(t3074);
t3079 = koppelP(6,3);
t3073 = xP(6);
t3037 = sin(t3073);
t3040 = cos(t3073);
t3085 = koppelP(6,1);
t3209 = koppelP(6,2);
t3220 = t3037 * t3209 - t3040 * t3085;
t2943 = t3038 * t3220 + t3041 * t3079;
t2982 = t3037 * t3085 + t3040 * t3209;
t3075 = xP(4);
t3039 = sin(t3075);
t3042 = cos(t3075);
t2886 = t2943 * t3039 - t2982 * t3042;
t2883 = t2943 * t3042 + t2982 * t3039;
t3080 = koppelP(5,3);
t3086 = koppelP(5,1);
t3210 = koppelP(5,2);
t3219 = t3037 * t3210 - t3040 * t3086;
t2945 = t3038 * t3219 + t3041 * t3080;
t2983 = t3037 * t3086 + t3040 * t3210;
t2887 = t2945 * t3039 - t2983 * t3042;
t2884 = t2945 * t3042 + t2983 * t3039;
t3081 = koppelP(4,3);
t3087 = koppelP(4,1);
t3211 = koppelP(4,2);
t3218 = t3037 * t3211 - t3040 * t3087;
t2947 = t3038 * t3218 + t3041 * t3081;
t2984 = t3037 * t3087 + t3040 * t3211;
t2888 = t2947 * t3039 - t2984 * t3042;
t2885 = t2947 * t3042 + t2984 * t3039;
t3082 = koppelP(3,3);
t3088 = koppelP(3,1);
t3212 = koppelP(3,2);
t3217 = t3037 * t3212 - t3040 * t3088;
t2949 = t3038 * t3217 + t3041 * t3082;
t2985 = t3037 * t3088 + t3040 * t3212;
t2895 = t2949 * t3039 - t2985 * t3042;
t2896 = t2949 * t3042 + t2985 * t3039;
t3083 = koppelP(2,3);
t3089 = koppelP(2,1);
t3213 = koppelP(2,2);
t3216 = t3037 * t3213 - t3040 * t3089;
t2951 = t3038 * t3216 + t3041 * t3083;
t2986 = t3037 * t3089 + t3040 * t3213;
t2899 = t2951 * t3039 - t2986 * t3042;
t2900 = t2951 * t3042 + t2986 * t3039;
t3084 = koppelP(1,3);
t3090 = koppelP(1,1);
t3214 = koppelP(1,2);
t3215 = t3037 * t3214 - t3040 * t3090;
t2953 = t3038 * t3215 + t3041 * t3084;
t2987 = t3037 * t3090 + t3040 * t3214;
t2903 = t2953 * t3039 - t2987 * t3042;
t2904 = t2953 * t3042 + t2987 * t3039;
t3076 = rSges(4,3);
t3077 = rSges(4,2);
t3078 = rSges(4,1);
t3097 = t3037 * t3077 - t3040 * t3078;
t3221 = -t3038 * t3076 + t3041 * t3097;
t3072 = m(2) + m(3);
t3066 = legFrame(6,2);
t3031 = cos(t3066);
t3208 = g(1) * t3031;
t3067 = legFrame(5,2);
t3032 = cos(t3067);
t3207 = g(1) * t3032;
t3068 = legFrame(4,2);
t3033 = cos(t3068);
t3206 = g(1) * t3033;
t3069 = legFrame(3,2);
t3034 = cos(t3069);
t3205 = g(1) * t3034;
t3070 = legFrame(2,2);
t3035 = cos(t3070);
t3204 = g(1) * t3035;
t3071 = legFrame(1,2);
t3036 = cos(t3071);
t3203 = g(1) * t3036;
t3202 = cos(pkin(8));
t3054 = legFrame(6,3);
t3001 = sin(t3054);
t3013 = cos(t3054);
t3048 = sin(pkin(8));
t2963 = t3001 * t3048 - t3013 * t3202;
t3049 = sin(pkin(5));
t3047 = sin(pkin(9));
t3051 = cos(pkin(9));
t3053 = cos(pkin(4));
t3143 = t3051 * t3053;
t2961 = t3047 * t3202 + t3048 * t3143;
t2962 = t3047 * t3048 - t3143 * t3202;
t3108 = t2961 * t3001 + t2962 * t3013;
t3050 = sin(pkin(4));
t3052 = cos(pkin(5));
t3144 = t3050 * t3052;
t3201 = (t2963 * t3144 + t3049 * t3108) * t3031;
t3055 = legFrame(5,3);
t3002 = sin(t3055);
t3014 = cos(t3055);
t2964 = t3002 * t3048 - t3014 * t3202;
t3106 = t2961 * t3002 + t2962 * t3014;
t3200 = (t2964 * t3144 + t3049 * t3106) * t3032;
t3056 = legFrame(4,3);
t3003 = sin(t3056);
t3015 = cos(t3056);
t2965 = t3003 * t3048 - t3015 * t3202;
t3104 = t2961 * t3003 + t2962 * t3015;
t3199 = (t2965 * t3144 + t3049 * t3104) * t3033;
t3057 = legFrame(3,3);
t3004 = sin(t3057);
t3016 = cos(t3057);
t2966 = t3004 * t3048 - t3016 * t3202;
t3102 = t2961 * t3004 + t2962 * t3016;
t3198 = (t2966 * t3144 + t3049 * t3102) * t3034;
t3058 = legFrame(2,3);
t3005 = sin(t3058);
t3017 = cos(t3058);
t2967 = t3005 * t3048 - t3017 * t3202;
t3100 = t2961 * t3005 + t2962 * t3017;
t3197 = (t2967 * t3144 + t3049 * t3100) * t3035;
t3059 = legFrame(1,3);
t3006 = sin(t3059);
t3018 = cos(t3059);
t2968 = t3006 * t3048 - t3018 * t3202;
t3098 = t2961 * t3006 + t2962 * t3018;
t3196 = (t2968 * t3144 + t3049 * t3098) * t3036;
t3195 = t2963 * t3031;
t3194 = t2963 * t3053;
t3193 = t2964 * t3032;
t3192 = t2964 * t3053;
t3191 = t2965 * t3033;
t3190 = t2965 * t3053;
t3189 = t2966 * t3034;
t3188 = t2966 * t3053;
t3187 = t2967 * t3035;
t3186 = t2967 * t3053;
t3185 = t2968 * t3036;
t3184 = t2968 * t3053;
t3060 = legFrame(6,1);
t3007 = sin(t3060);
t3025 = sin(t3066);
t3171 = t3007 * t3025;
t3061 = legFrame(5,1);
t3008 = sin(t3061);
t3026 = sin(t3067);
t3170 = t3008 * t3026;
t3062 = legFrame(4,1);
t3009 = sin(t3062);
t3027 = sin(t3068);
t3169 = t3009 * t3027;
t3063 = legFrame(3,1);
t3010 = sin(t3063);
t3028 = sin(t3069);
t3168 = t3010 * t3028;
t3064 = legFrame(2,1);
t3011 = sin(t3064);
t3029 = sin(t3070);
t3167 = t3011 * t3029;
t3065 = legFrame(1,1);
t3012 = sin(t3065);
t3030 = sin(t3071);
t3166 = t3012 * t3030;
t3019 = cos(t3060);
t3165 = t3019 * t3025;
t3020 = cos(t3061);
t3164 = t3020 * t3026;
t3021 = cos(t3062);
t3163 = t3021 * t3027;
t3022 = cos(t3063);
t3162 = t3022 * t3028;
t3023 = cos(t3064);
t3161 = t3023 * t3029;
t3024 = cos(t3065);
t3160 = t3024 * t3030;
t3158 = t3038 * t3077;
t3157 = t3038 * t3078;
t3150 = t3041 * t3076;
t3043 = 0.1e1 / t3047;
t3044 = 0.1e1 / t3049;
t3149 = t3043 * t3044;
t3045 = 0.1e1 / t3050;
t3148 = t3045 * (m(1) + t3072);
t3147 = t3045 * t3072;
t3146 = t3047 * t3050;
t3145 = t3047 * t3053;
t3136 = m(3) * t3149;
t3135 = t3043 * t3148;
t2865 = -t3001 * t3208 + (-t3001 * t3171 + t3013 * t3019) * g(2) + (t3001 * t3165 + t3007 * t3013) * g(3);
t2866 = t3013 * t3208 + (t3001 * t3019 + t3013 * t3171) * g(2) + (t3001 * t3007 - t3013 * t3165) * g(3);
t2841 = t2865 * t3202 - t3048 * t2866;
t2937 = g(1) * t3025 + (-g(2) * t3007 + g(3) * t3019) * t3031;
t2823 = t2841 * t3050 - t2937 * t3053;
t2817 = ((t2841 * t3053 + t2937 * t3050) * t3051 + (-t2865 * t3048 - t2866 * t3202) * t3047) * t3049 + t3052 * t2823;
t3134 = t2817 * t3136;
t2867 = -t3002 * t3207 + (-t3002 * t3170 + t3014 * t3020) * g(2) + (t3002 * t3164 + t3008 * t3014) * g(3);
t2868 = t3014 * t3207 + (t3002 * t3020 + t3014 * t3170) * g(2) + (t3002 * t3008 - t3014 * t3164) * g(3);
t2842 = t2867 * t3202 - t3048 * t2868;
t2938 = g(1) * t3026 + (-g(2) * t3008 + g(3) * t3020) * t3032;
t2824 = t2842 * t3050 - t2938 * t3053;
t2818 = ((t2842 * t3053 + t2938 * t3050) * t3051 + (-t2867 * t3048 - t2868 * t3202) * t3047) * t3049 + t3052 * t2824;
t3133 = t2818 * t3136;
t2869 = -t3003 * t3206 + (-t3003 * t3169 + t3015 * t3021) * g(2) + (t3003 * t3163 + t3009 * t3015) * g(3);
t2870 = t3015 * t3206 + (t3003 * t3021 + t3015 * t3169) * g(2) + (t3003 * t3009 - t3015 * t3163) * g(3);
t2843 = t2869 * t3202 - t3048 * t2870;
t2939 = g(1) * t3027 + (-g(2) * t3009 + g(3) * t3021) * t3033;
t2825 = t2843 * t3050 - t2939 * t3053;
t2819 = ((t2843 * t3053 + t2939 * t3050) * t3051 + (-t2869 * t3048 - t2870 * t3202) * t3047) * t3049 + t3052 * t2825;
t3132 = t2819 * t3136;
t2871 = -t3004 * t3205 + (-t3004 * t3168 + t3016 * t3022) * g(2) + (t3004 * t3162 + t3010 * t3016) * g(3);
t2872 = t3016 * t3205 + (t3004 * t3022 + t3016 * t3168) * g(2) + (t3004 * t3010 - t3016 * t3162) * g(3);
t2844 = t2871 * t3202 - t3048 * t2872;
t2940 = g(1) * t3028 + (-g(2) * t3010 + g(3) * t3022) * t3034;
t2826 = t2844 * t3050 - t2940 * t3053;
t2820 = ((t2844 * t3053 + t2940 * t3050) * t3051 + (-t2871 * t3048 - t2872 * t3202) * t3047) * t3049 + t3052 * t2826;
t3131 = t2820 * t3136;
t2873 = -t3005 * t3204 + (-t3005 * t3167 + t3017 * t3023) * g(2) + (t3005 * t3161 + t3011 * t3017) * g(3);
t2874 = t3017 * t3204 + (t3005 * t3023 + t3017 * t3167) * g(2) + (t3005 * t3011 - t3017 * t3161) * g(3);
t2845 = t2873 * t3202 - t3048 * t2874;
t2941 = g(1) * t3029 + (-g(2) * t3011 + g(3) * t3023) * t3035;
t2827 = t2845 * t3050 - t2941 * t3053;
t2821 = ((t2845 * t3053 + t2941 * t3050) * t3051 + (-t2873 * t3048 - t2874 * t3202) * t3047) * t3049 + t3052 * t2827;
t3130 = t2821 * t3136;
t2875 = -t3006 * t3203 + (-t3006 * t3166 + t3018 * t3024) * g(2) + (t3006 * t3160 + t3012 * t3018) * g(3);
t2876 = t3018 * t3203 + (t3006 * t3024 + t3018 * t3166) * g(2) + (t3006 * t3012 - t3018 * t3160) * g(3);
t2846 = t2875 * t3202 - t3048 * t2876;
t2942 = g(1) * t3030 + (-g(2) * t3012 + g(3) * t3024) * t3036;
t2828 = t2846 * t3050 - t2942 * t3053;
t2822 = ((t2846 * t3053 + t2942 * t3050) * t3051 + (-t2875 * t3048 - t2876 * t3202) * t3047) * t3049 + t3052 * t2828;
t3129 = t2822 * t3136;
t3128 = t2937 * t3135;
t3127 = t2938 * t3135;
t3126 = t2939 * t3135;
t3125 = t2940 * t3135;
t3124 = t2941 * t3135;
t3123 = t2942 * t3135;
t3122 = t3147 * t3149;
t2969 = t3001 * t3202 + t3013 * t3048;
t3121 = t2963 * t3165 + t3007 * t2969;
t2970 = t3002 * t3202 + t3014 * t3048;
t3120 = t2964 * t3164 + t3008 * t2970;
t2971 = t3003 * t3202 + t3015 * t3048;
t3119 = t2965 * t3163 + t3009 * t2971;
t2972 = t3004 * t3202 + t3016 * t3048;
t3118 = t2966 * t3162 + t3010 * t2972;
t2973 = t3005 * t3202 + t3017 * t3048;
t3117 = t2967 * t3161 + t3011 * t2973;
t2974 = t3006 * t3202 + t3018 * t3048;
t3116 = t2968 * t3160 + t3012 * t2974;
t3115 = t2823 * t3122;
t3114 = t2824 * t3122;
t3113 = t2825 * t3122;
t3112 = t2826 * t3122;
t3111 = t2827 * t3122;
t3110 = t2828 * t3122;
t3109 = t2961 * t3013 - t2962 * t3001;
t3107 = t2961 * t3014 - t2962 * t3002;
t3105 = t2961 * t3015 - t2962 * t3003;
t3103 = t2961 * t3016 - t2962 * t3004;
t3101 = t2961 * t3017 - t2962 * t3005;
t3099 = t2961 * t3018 - t2962 * t3006;
t3096 = t2963 * t3171 - t2969 * t3019;
t3095 = t2964 * t3170 - t2970 * t3020;
t3094 = t2965 * t3169 - t2971 * t3021;
t3093 = t2966 * t3168 - t2972 * t3022;
t3092 = t2967 * t3167 - t2973 * t3023;
t3091 = t2968 * t3166 - t2974 * t3024;
t2936 = t3038 * t3084 - t3041 * t3215;
t2935 = t3038 * t3083 - t3041 * t3216;
t2934 = t3038 * t3082 - t3041 * t3217;
t2933 = t3038 * t3081 - t3041 * t3218;
t2932 = t3038 * t3080 - t3041 * t3219;
t2931 = t3038 * t3079 - t3041 * t3220;
t2930 = t2974 * t3030 * t3053 + t3050 * t3036;
t2929 = t2973 * t3029 * t3053 + t3050 * t3035;
t2928 = t2972 * t3028 * t3053 + t3050 * t3034;
t2927 = t2971 * t3027 * t3053 + t3050 * t3033;
t2926 = t2970 * t3026 * t3053 + t3050 * t3032;
t2925 = t2969 * t3025 * t3053 + t3050 * t3031;
t2882 = (-t2968 * t3051 - t2974 * t3145) * t3036 + t3030 * t3146;
t2881 = (-t2967 * t3051 - t2973 * t3145) * t3035 + t3029 * t3146;
t2880 = (-t2966 * t3051 - t2972 * t3145) * t3034 + t3028 * t3146;
t2879 = (-t2965 * t3051 - t2971 * t3145) * t3033 + t3027 * t3146;
t2878 = (-t2964 * t3051 - t2970 * t3145) * t3032 + t3026 * t3146;
t2877 = (-t2963 * t3051 - t2969 * t3145) * t3031 + t3025 * t3146;
t2858 = (-t2930 * t3012 - t3024 * t3184) * t3047 - t3091 * t3051;
t2857 = (t2930 * t3024 - t3012 * t3184) * t3047 + t3051 * t3116;
t2856 = (-t2929 * t3011 - t3023 * t3186) * t3047 - t3092 * t3051;
t2855 = (t2929 * t3023 - t3011 * t3186) * t3047 + t3051 * t3117;
t2854 = (-t2928 * t3010 - t3022 * t3188) * t3047 - t3093 * t3051;
t2853 = (t2928 * t3022 - t3010 * t3188) * t3047 + t3051 * t3118;
t2852 = (-t2927 * t3009 - t3021 * t3190) * t3047 - t3094 * t3051;
t2851 = (t2927 * t3021 - t3009 * t3190) * t3047 + t3051 * t3119;
t2850 = (-t2926 * t3008 - t3020 * t3192) * t3047 - t3095 * t3051;
t2849 = (t2926 * t3020 - t3008 * t3192) * t3047 + t3051 * t3120;
t2848 = (-t2925 * t3007 - t3019 * t3194) * t3047 - t3096 * t3051;
t2847 = (t2925 * t3019 - t3007 * t3194) * t3047 + t3051 * t3121;
t2840 = (-t3024 * t3099 + t3098 * t3166) * t3049 + t3091 * t3144;
t2839 = (-t3023 * t3101 + t3100 * t3167) * t3049 + t3092 * t3144;
t2838 = (-t3022 * t3103 + t3102 * t3168) * t3049 + t3093 * t3144;
t2837 = (-t3021 * t3105 + t3104 * t3169) * t3049 + t3094 * t3144;
t2836 = (-t3020 * t3107 + t3106 * t3170) * t3049 + t3095 * t3144;
t2835 = (-t3019 * t3109 + t3108 * t3171) * t3049 + t3096 * t3144;
t2834 = (-t3012 * t3099 - t3098 * t3160) * t3049 - t3116 * t3144;
t2833 = (-t3011 * t3101 - t3100 * t3161) * t3049 - t3117 * t3144;
t2832 = (-t3010 * t3103 - t3102 * t3162) * t3049 - t3118 * t3144;
t2831 = (-t3009 * t3105 - t3104 * t3163) * t3049 - t3119 * t3144;
t2830 = (-t3008 * t3107 - t3106 * t3164) * t3049 - t3120 * t3144;
t2829 = (-t3007 * t3109 - t3108 * t3165) * t3049 - t3121 * t3144;
t1 = [-m(4) * g(1) + ((-t2877 * t2937 - t2878 * t2938 - t2879 * t2939 - t2880 * t2940 - t2881 * t2941 - t2882 * t2942) * t3148 + ((t2823 * t3201 + t2824 * t3200 + t2825 * t3199 + t2826 * t3198 + t2827 * t3197 + t2828 * t3196) * t3147 + (-t2817 * t3195 - t2818 * t3193 - t2819 * t3191 - t2820 * t3189 - t2821 * t3187 - t2822 * t3185) * m(3)) * t3044) * t3043; -m(4) * g(2) + ((-t2848 * t2937 - t2850 * t2938 - t2852 * t2939 - t2854 * t2940 - t2856 * t2941 - t2858 * t2942) * t3148 + ((t2823 * t2835 + t2824 * t2836 + t2825 * t2837 + t2826 * t2838 + t2827 * t2839 + t2828 * t2840) * t3147 + (-t2817 * t3096 - t2818 * t3095 - t2819 * t3094 - t2820 * t3093 - t2821 * t3092 - t2822 * t3091) * m(3)) * t3044) * t3043; -m(4) * g(3) + ((-t2847 * t2937 - t2849 * t2938 - t2851 * t2939 - t2853 * t2940 - t2855 * t2941 - t2857 * t2942) * t3148 + ((t2823 * t2829 + t2824 * t2830 + t2825 * t2831 + t2826 * t2832 + t2827 * t2833 + t2828 * t2834) * t3147 + (t2817 * t3121 + t2818 * t3120 + t2819 * t3119 + t2820 * t3118 + t2821 * t3117 + t2822 * t3116) * m(3)) * t3044) * t3043; -(-t2857 * t2903 - t2858 * t2904) * t3123 + (-t2834 * t2903 - t2840 * t2904) * t3110 + (-t2903 * t3116 + t2904 * t3091) * t3129 - (-t2855 * t2899 - t2856 * t2900) * t3124 + (-t2833 * t2899 - t2839 * t2900) * t3111 + (-t2899 * t3117 + t2900 * t3092) * t3130 - (-t2853 * t2895 - t2854 * t2896) * t3125 + (-t2832 * t2895 - t2838 * t2896) * t3112 + (-t2895 * t3118 + t2896 * t3093) * t3131 - (-t2851 * t2888 - t2852 * t2885) * t3126 + (-t2831 * t2888 - t2837 * t2885) * t3113 + (t2885 * t3094 - t2888 * t3119) * t3132 - (-t2849 * t2887 - t2850 * t2884) * t3127 + (-t2830 * t2887 - t2836 * t2884) * t3114 + (t2884 * t3095 - t2887 * t3120) * t3133 - (-t2847 * t2886 - t2848 * t2883) * t3128 + (-t2829 * t2886 - t2835 * t2883) * t3115 + (t2883 * t3096 - t2886 * t3121) * t3134 + (((-g(2) * t3157 - g(3) * t3077) * t3040 + g(2) * t3150) * t3042 + t3039 * ((g(2) * t3077 - g(3) * t3157) * t3040 + g(3) * t3150) + ((g(2) * t3158 - g(3) * t3078) * t3042 + t3039 * (g(2) * t3078 + g(3) * t3158)) * t3037) * m(4); -(-t2857 * t2936 + t2882 * t2904) * t3123 + (-t2834 * t2936 + t2904 * t3196) * t3110 + (-t2904 * t3185 - t2936 * t3116) * t3129 - (-t2855 * t2935 + t2881 * t2900) * t3124 + (-t2833 * t2935 + t2900 * t3197) * t3111 + (-t2900 * t3187 - t2935 * t3117) * t3130 - (-t2853 * t2934 + t2880 * t2896) * t3125 + (-t2832 * t2934 + t2896 * t3198) * t3112 + (-t2896 * t3189 - t2934 * t3118) * t3131 - (-t2851 * t2933 + t2879 * t2885) * t3126 + (-t2831 * t2933 + t2885 * t3199) * t3113 + (-t2885 * t3191 - t2933 * t3119) * t3132 - (-t2849 * t2932 + t2878 * t2884) * t3127 + (-t2830 * t2932 + t2884 * t3200) * t3114 + (-t2884 * t3193 - t2932 * t3120) * t3133 - (-t2847 * t2931 + t2877 * t2883) * t3128 + (-t2829 * t2931 + t2883 * t3201) * t3115 + (-t2883 * t3195 - t2931 * t3121) * t3134 - (t3221 * g(3) + ((t3038 * t3097 + t3150) * t3042 + t3039 * (t3037 * t3078 + t3040 * t3077)) * g(1)) * m(4); -(t2858 * t2936 + t2882 * t2903) * t3123 + (t2840 * t2936 + t2903 * t3196) * t3110 + (-t2903 * t3185 - t2936 * t3091) * t3129 - (t2856 * t2935 + t2881 * t2899) * t3124 + (t2839 * t2935 + t2899 * t3197) * t3111 + (-t2899 * t3187 - t2935 * t3092) * t3130 - (t2854 * t2934 + t2880 * t2895) * t3125 + (t2838 * t2934 + t2895 * t3198) * t3112 + (-t2895 * t3189 - t2934 * t3093) * t3131 - (t2852 * t2933 + t2879 * t2888) * t3126 + (t2837 * t2933 + t2888 * t3199) * t3113 + (-t2888 * t3191 - t2933 * t3094) * t3132 - (t2850 * t2932 + t2878 * t2887) * t3127 + (t2836 * t2932 + t2887 * t3200) * t3114 + (-t2887 * t3193 - t2932 * t3095) * t3133 - (t2848 * t2931 + t2877 * t2886) * t3128 + (t2835 * t2931 + t2886 * t3201) * t3115 + (-t2886 * t3195 - t2931 * t3096) * t3134 + (t3221 * g(2) + (-t3039 * t3150 + (t3039 * t3157 + t3042 * t3077) * t3040 + (-t3039 * t3158 + t3042 * t3078) * t3037) * g(1)) * m(4);];
taugX  = t1;
