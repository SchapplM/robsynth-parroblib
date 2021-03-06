% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR9V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRRR9V1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:52:44
% EndTime: 2020-08-06 18:52:57
% DurationCPUTime: 12.33s
% Computational Cost: add. (11727->564), mult. (17498->1078), div. (2064->17), fcn. (15720->26), ass. (0->448)
t2903 = sin(pkin(7));
t2904 = cos(pkin(7));
t3230 = MDP(4) * t2904 - MDP(5) * t2903;
t3223 = 2 * pkin(3);
t3222 = 2 * MDP(6);
t3221 = 4 * MDP(9);
t2893 = t2904 ^ 2;
t3220 = 0.2e1 * t2893;
t3219 = -0.2e1 * t2903;
t2907 = (pkin(5) + qJ(2,1));
t2906 = (pkin(5) + qJ(2,2));
t2905 = (pkin(5) + qJ(2,3));
t2917 = cos(qJ(3,3));
t2897 = t2917 ^ 2;
t3218 = t2897 * pkin(3);
t2919 = cos(qJ(3,2));
t2899 = t2919 ^ 2;
t3217 = t2899 * pkin(3);
t2921 = cos(qJ(3,1));
t2901 = t2921 ^ 2;
t3216 = t2901 * pkin(3);
t3215 = t2917 * pkin(2);
t3214 = t2919 * pkin(2);
t3213 = t2921 * pkin(2);
t2925 = 1 / pkin(3);
t3210 = MDP(10) * t2925;
t3209 = MDP(11) * t2925;
t3208 = MDP(12) / pkin(3) ^ 2;
t2912 = sin(qJ(1,3));
t2890 = pkin(6) + t2905;
t2918 = cos(qJ(1,3));
t3047 = pkin(1) * t2912 - t2890 * t2918;
t2911 = sin(qJ(3,3));
t3125 = t2903 * t2911;
t2807 = t3047 * t3125 + (t2897 - 0.1e1) * t2912 * pkin(3);
t2816 = pkin(1) * t2911 + (-pkin(3) + t3215 + 0.2e1 * t3218) * t2903;
t2923 = -pkin(3) / 0.2e1;
t2846 = t3218 + t3215 / 0.2e1 + t2923;
t2908 = legFrame(3,2);
t2878 = sin(t2908);
t2881 = cos(t2908);
t3057 = t2912 * t3125;
t2959 = pkin(2) * t3057 + (t3057 * t3223 - t3047) * t2917;
t3132 = t2878 * t2912;
t2924 = pkin(2) / 0.2e1;
t3159 = (pkin(3) * t2917 + t2924) * t2911;
t2871 = pkin(1) * t2903;
t3162 = (-pkin(3) * t2911 + t2871) * t2917;
t2763 = (-t2846 * t3132 + t2881 * t3159) * t3220 + (t2881 * t2816 + t2959 * t2878) * t2904 + t2807 * t2878 + t2881 * t3162;
t3119 = t2904 * t2917;
t2841 = t3119 - t3125;
t2835 = 0.1e1 / t2841;
t3207 = t2763 * t2835;
t2873 = 0.1e1 / t2890 ^ 2;
t3206 = t2763 * t2873;
t3129 = t2881 * t2912;
t2764 = (t2846 * t3129 + t2878 * t3159) * t3220 + (t2878 * t2816 - t2959 * t2881) * t2904 - t2807 * t2881 + t2878 * t3162;
t3205 = t2764 * t2835;
t3204 = t2764 * t2873;
t2914 = sin(qJ(1,2));
t2891 = pkin(6) + t2906;
t2920 = cos(qJ(1,2));
t3046 = pkin(1) * t2914 - t2891 * t2920;
t2913 = sin(qJ(3,2));
t3124 = t2903 * t2913;
t2808 = t3046 * t3124 + (t2899 - 0.1e1) * t2914 * pkin(3);
t2817 = pkin(1) * t2913 + (-pkin(3) + t3214 + 0.2e1 * t3217) * t2903;
t2847 = t3217 + t3214 / 0.2e1 + t2923;
t2909 = legFrame(2,2);
t2879 = sin(t2909);
t2882 = cos(t2909);
t3056 = t2914 * t3124;
t2958 = pkin(2) * t3056 + (t3056 * t3223 - t3046) * t2919;
t3131 = t2879 * t2914;
t3158 = (pkin(3) * t2919 + t2924) * t2913;
t3161 = (-pkin(3) * t2913 + t2871) * t2919;
t2765 = (-t2847 * t3131 + t2882 * t3158) * t3220 + (t2882 * t2817 + t2958 * t2879) * t2904 + t2808 * t2879 + t2882 * t3161;
t3118 = t2904 * t2919;
t2842 = t3118 - t3124;
t2836 = 0.1e1 / t2842;
t3203 = t2765 * t2836;
t2875 = 0.1e1 / t2891 ^ 2;
t3202 = t2765 * t2875;
t3128 = t2882 * t2914;
t2766 = (t2847 * t3128 + t2879 * t3158) * t3220 + (t2879 * t2817 - t2958 * t2882) * t2904 - t2808 * t2882 + t2879 * t3161;
t3201 = t2766 * t2836;
t3200 = t2766 * t2875;
t2916 = sin(qJ(1,1));
t2892 = pkin(6) + t2907;
t2922 = cos(qJ(1,1));
t3048 = t2916 * pkin(1) - t2892 * t2922;
t2915 = sin(qJ(3,1));
t3123 = t2903 * t2915;
t2809 = t3048 * t3123 + (t2901 - 0.1e1) * t2916 * pkin(3);
t2818 = pkin(1) * t2915 + (-pkin(3) + t3213 + 0.2e1 * t3216) * t2903;
t2848 = t3216 + t3213 / 0.2e1 + t2923;
t2910 = legFrame(1,2);
t2880 = sin(t2910);
t2883 = cos(t2910);
t3055 = t2916 * t3123;
t2957 = pkin(2) * t3055 + (t3055 * t3223 - t3048) * t2921;
t3130 = t2880 * t2916;
t3157 = (pkin(3) * t2921 + t2924) * t2915;
t3160 = (-pkin(3) * t2915 + t2871) * t2921;
t2767 = (-t2848 * t3130 + t2883 * t3157) * t3220 + (t2883 * t2818 + t2957 * t2880) * t2904 + t2809 * t2880 + t2883 * t3160;
t3117 = t2904 * t2921;
t2840 = t3117 - t3123;
t2834 = 0.1e1 / t2840;
t3199 = t2767 * t2834;
t2877 = 0.1e1 / t2892 ^ 2;
t3198 = t2767 * t2877;
t3127 = t2883 * t2916;
t2768 = (t2848 * t3127 + t2880 * t3157) * t3220 + (t2880 * t2818 - t2957 * t2883) * t2904 - t2809 * t2883 + t2880 * t3160;
t3197 = t2768 * t2834;
t3196 = t2768 * t2877;
t2858 = pkin(2) * t2904 + pkin(1);
t2894 = pkin(7) + qJ(3,3);
t2868 = cos(t2894);
t2813 = t2890 * t2912 + (pkin(3) * t2868 + t2858) * t2918;
t2872 = 0.1e1 / t2890;
t2804 = t2813 * t2872;
t3141 = t2872 * t2918;
t3102 = pkin(1) * t3141;
t2793 = t2804 - t3102;
t3195 = t2793 * t2835;
t2895 = pkin(7) + qJ(3,2);
t2869 = cos(t2895);
t2814 = t2891 * t2914 + (pkin(3) * t2869 + t2858) * t2920;
t2874 = 0.1e1 / t2891;
t2805 = t2814 * t2874;
t3138 = t2874 * t2920;
t3101 = pkin(1) * t3138;
t2795 = t2805 - t3101;
t3194 = t2795 * t2836;
t2896 = pkin(7) + qJ(3,1);
t2870 = cos(t2896);
t2815 = t2892 * t2916 + (pkin(3) * t2870 + t2858) * t2922;
t2876 = 0.1e1 / t2892;
t2806 = t2815 * t2876;
t3135 = t2876 * t2922;
t3100 = pkin(1) * t3135;
t2797 = t2806 - t3100;
t3193 = t2797 * t2834;
t3192 = t2813 * t2873;
t3191 = t2814 * t2875;
t3190 = t2815 * t2877;
t2865 = sin(t2894);
t2825 = t2865 * t2881 - t2868 * t3132;
t2859 = 0.1e1 / t2868;
t3189 = t2825 * t2859;
t3188 = t2825 * t2872;
t3187 = t2825 * t2878;
t2826 = t2865 * t2878 + t2868 * t3129;
t3186 = t2826 * t2859;
t3185 = t2826 * t2872;
t3184 = t2826 * t2881;
t2866 = sin(t2895);
t2827 = t2866 * t2882 - t2869 * t3131;
t2861 = 0.1e1 / t2869;
t3183 = t2827 * t2861;
t3182 = t2827 * t2874;
t3181 = t2827 * t2879;
t2828 = t2866 * t2879 + t2869 * t3128;
t3180 = t2828 * t2861;
t3179 = t2828 * t2874;
t3178 = t2828 * t2882;
t2867 = sin(t2896);
t2829 = t2867 * t2883 - t2870 * t3130;
t2863 = 0.1e1 / t2870;
t3177 = t2829 * t2863;
t3176 = t2829 * t2876;
t3175 = t2829 * t2880;
t2830 = t2867 * t2880 + t2870 * t3127;
t3174 = t2830 * t2863;
t3173 = t2830 * t2876;
t3172 = t2830 * t2883;
t3171 = t2834 * t2876;
t3170 = t2834 * t2877;
t3169 = t2835 * t2872;
t3168 = t2835 * t2873;
t3167 = t2836 * t2874;
t3166 = t2836 * t2875;
t3122 = t2903 * t2917;
t2843 = t2904 * t2911 + t3122;
t3165 = t2843 * t2873;
t3121 = t2903 * t2919;
t2844 = t2904 * t2913 + t3121;
t3164 = t2844 * t2875;
t3120 = t2903 * t2921;
t2845 = t2904 * t2915 + t3120;
t3163 = t2845 * t2877;
t3156 = t2859 * t2872;
t2927 = pkin(1) ^ 2;
t2884 = qJ(2,3) ^ 2 + t2927;
t3155 = t2859 * t2884;
t2860 = 0.1e1 / t2868 ^ 2;
t3154 = t2860 * t2872;
t3153 = t2860 * t2873;
t3152 = t2860 * t2881;
t3151 = t2861 * t2874;
t2885 = qJ(2,2) ^ 2 + t2927;
t3150 = t2861 * t2885;
t2862 = 0.1e1 / t2869 ^ 2;
t3149 = t2862 * t2874;
t3148 = t2862 * t2875;
t3147 = t2862 * t2882;
t3146 = t2863 * t2876;
t2886 = qJ(2,1) ^ 2 + t2927;
t3145 = t2863 * t2886;
t2864 = 0.1e1 / t2870 ^ 2;
t3144 = t2864 * t2876;
t3143 = t2864 * t2877;
t3142 = t2864 * t2883;
t3140 = t2873 * t2918 ^ 2;
t3139 = t2873 * t2918;
t3137 = t2875 * t2920 ^ 2;
t3136 = t2875 * t2920;
t3134 = t2877 * t2922 ^ 2;
t3133 = t2877 * t2922;
t3126 = t2903 * t2904;
t3116 = t2905 * t2925;
t3115 = t2906 * t2925;
t3114 = t2907 * t2925;
t3113 = t2893 - 0.1e1 / 0.2e1;
t3112 = pkin(2) * t3220;
t3111 = pkin(2) * t3219;
t3110 = 2 * t3210;
t3109 = 2 * t3209;
t3108 = -0.2e1 * t3125;
t3107 = -0.2e1 * t3124;
t3106 = -0.2e1 * t3123;
t3105 = pkin(2) * t3141;
t3104 = pkin(2) * t3138;
t3103 = pkin(2) * t3135;
t3099 = t2763 * t3168;
t3098 = t2764 * t3168;
t3097 = t2765 * t3166;
t3096 = t2766 * t3166;
t3095 = t2767 * t3170;
t3094 = t2768 * t3170;
t3093 = t2859 * t3192;
t3092 = t2861 * t3191;
t3091 = t2863 * t3190;
t2819 = t2825 ^ 2;
t3090 = t2819 * t3153;
t2820 = t2826 ^ 2;
t3089 = t2820 * t3153;
t2821 = t2827 ^ 2;
t3088 = t2821 * t3148;
t2822 = t2828 ^ 2;
t3087 = t2822 * t3148;
t2823 = t2829 ^ 2;
t3086 = t2823 * t3143;
t2824 = t2830 ^ 2;
t3085 = t2824 * t3143;
t3084 = t2825 * t3156;
t3083 = t2826 * t3156;
t3082 = t2827 * t3151;
t3081 = t2828 * t3151;
t3080 = t2829 * t3146;
t3079 = t2830 * t3146;
t3078 = t2834 * t3163;
t3077 = t2834 * t3133;
t3076 = t2835 * t3165;
t3075 = t2835 * t3139;
t3074 = t2836 * t3164;
t3073 = t2836 * t3136;
t2837 = t2843 ^ 2;
t3072 = t2837 * t3153;
t2838 = t2844 ^ 2;
t3071 = t2838 * t3148;
t2839 = t2845 ^ 2;
t3070 = t2839 * t3143;
t3069 = t2859 * t3141;
t3068 = t2859 * t3139;
t3067 = t2905 * t3154;
t3066 = t2861 * t3138;
t3065 = t2861 * t3136;
t3064 = t2906 * t3149;
t3063 = t2863 * t3135;
t3062 = t2863 * t3133;
t3061 = t2907 * t3144;
t3060 = t2903 * t3116;
t3059 = t2903 * t3115;
t3058 = t2903 * t3114;
t3054 = t2872 * t3112;
t3053 = t2872 * t3111;
t3052 = t2874 * t3112;
t3051 = t2874 * t3111;
t3050 = t2876 * t3112;
t3049 = t2876 * t3111;
t3045 = t2911 * t3105;
t3044 = t2917 * t3105;
t3043 = t2913 * t3104;
t3042 = t2919 * t3104;
t3041 = t2915 * t3103;
t3040 = t2921 * t3103;
t3039 = pkin(1) * t3084;
t3038 = pkin(1) * t3083;
t3037 = pkin(1) * t3082;
t3036 = pkin(1) * t3081;
t3035 = pkin(1) * t3080;
t3034 = pkin(1) * t3079;
t3033 = t2763 * t3076;
t3032 = t2764 * t3076;
t3031 = t2765 * t3074;
t3030 = t2766 * t3074;
t3029 = t2767 * t3078;
t3028 = t2768 * t3078;
t3027 = t2841 * t3093;
t3026 = t2843 * t3093;
t3025 = t2842 * t3092;
t3024 = t2844 * t3092;
t3023 = t2840 * t3091;
t3022 = t2845 * t3091;
t3021 = t2825 * t2826 * t3153;
t3020 = t2825 * t3068;
t3019 = t3152 * t3188;
t3018 = t2826 * t3068;
t3017 = t2826 * t2878 * t3154;
t3016 = t2827 * t2828 * t3148;
t3015 = t2827 * t3065;
t3014 = t3147 * t3182;
t3013 = t2828 * t3065;
t3012 = t2828 * t2879 * t3149;
t3011 = t2829 * t2830 * t3143;
t3010 = t2829 * t3062;
t3009 = t3142 * t3176;
t3008 = t2830 * t3062;
t3007 = t2830 * t2880 * t3144;
t3006 = t2845 * t3077;
t3005 = t2843 * t3075;
t3004 = t2844 * t3073;
t3003 = t2837 * t3068;
t3002 = t2838 * t3065;
t3001 = t2839 * t3062;
t3000 = t2878 * t3069;
t2999 = t2881 * t3069;
t2998 = t2905 * t3069;
t2997 = t2925 * t3067;
t2996 = t2879 * t3066;
t2995 = t2882 * t3066;
t2994 = t2906 * t3066;
t2993 = t2925 * t3064;
t2992 = t2880 * t3063;
t2991 = t2883 * t3063;
t2990 = t2907 * t3063;
t2989 = t2925 * t3061;
t2988 = t3067 * t3187;
t2987 = t3067 * t3184;
t2986 = t3064 * t3181;
t2985 = t3064 * t3178;
t2984 = t3061 * t3175;
t2983 = t3061 * t3172;
t2982 = t2878 * t2998;
t2981 = t2881 * t2998;
t2980 = t2878 * t2997;
t2979 = t2881 * t2997;
t2978 = t2879 * t2994;
t2977 = t2882 * t2994;
t2976 = t2879 * t2993;
t2975 = t2882 * t2993;
t2974 = t2880 * t2990;
t2973 = t2883 * t2990;
t2972 = t2880 * t2989;
t2971 = t2883 * t2989;
t2810 = (t2897 - 0.1e1 / 0.2e1) * t3126 + t3113 * t2917 * t2911;
t2811 = (t2899 - 0.1e1 / 0.2e1) * t3126 + t3113 * t2919 * t2913;
t2812 = (t2901 - 0.1e1 / 0.2e1) * t3126 + t3113 * t2921 * t2915;
t2967 = (t2810 * t3020 + t2811 * t3015 + t2812 * t3010) * t3221 + (t2825 * t3003 + t2827 * t3002 + t2829 * t3001) * MDP(8) + (qJ(2,1) * t3010 + qJ(2,2) * t3015 + qJ(2,3) * t3020) * t3222 + (t3010 + t3015 + t3020) * MDP(1) + (t2840 * t2991 + t2841 * t2999 + t2842 * t2995) * t3209 + (t2843 * t2999 + t2844 * t2995 + t2845 * t2991) * t3210;
t2966 = (t2810 * t3018 + t2811 * t3013 + t2812 * t3008) * t3221 + (t2826 * t3003 + t2828 * t3002 + t2830 * t3001) * MDP(8) + (qJ(2,1) * t3008 + qJ(2,2) * t3013 + qJ(2,3) * t3018) * t3222 + (t3008 + t3013 + t3018) * MDP(1) + (t2840 * t2992 + t2841 * t3000 + t2842 * t2996) * t3209 + (t2843 * t3000 + t2844 * t2996 + t2845 * t2992) * t3210;
t2956 = 0.2e1 * t2843;
t2955 = 0.2e1 * t2844;
t2954 = 0.2e1 * t2845;
t2950 = (t3172 + t3175) * t3144;
t2951 = (t3178 + t3181) * t3149;
t2952 = (t3184 + t3187) * t3154;
t2953 = (t2840 * t2950 + t2841 * t2952 + t2842 * t2951) * t3209 + (t2843 * t2952 + t2844 * t2951 + t2845 * t2950) * t3210 + (t2810 * t3021 + t2811 * t3016 + t2812 * t3011) * t3221 + (t2837 * t3021 + t2838 * t3016 + t2839 * t3011) * MDP(8) + (qJ(2,1) * t3011 + qJ(2,2) * t3016 + qJ(2,3) * t3021) * t3222 + (t3011 + t3016 + t3021) * MDP(1) + (t2878 * t3152 + t2879 * t3147 + t2880 * t3142) * t3208;
t2943 = t2825 * t3054 - t2881 * t3060;
t2942 = t2826 * t3054 - t2878 * t3060;
t2941 = t2827 * t3052 - t2882 * t3059;
t2940 = t2828 * t3052 - t2879 * t3059;
t2939 = t2829 * t3050 - t2883 * t3058;
t2938 = t2830 * t3050 - t2880 * t3058;
t2937 = t2904 * (t2825 * t3053 - t2881 * t3116);
t2936 = t2904 * (t2826 * t3053 - t2878 * t3116);
t2935 = t2904 * (t2827 * t3051 - t2882 * t3115);
t2934 = t2904 * (t2828 * t3051 - t2879 * t3115);
t2933 = t2904 * (t2829 * t3049 - t2883 * t3114);
t2932 = t2904 * (t2830 * t3049 - t2880 * t3114);
t2796 = -t2806 + 0.2e1 * t3100;
t2794 = -t2805 + 0.2e1 * t3101;
t2792 = -t2804 + 0.2e1 * t3102;
t2791 = t3100 - t2806 / 0.2e1;
t2790 = t3101 - t2805 / 0.2e1;
t2789 = t3102 - t2804 / 0.2e1;
t2788 = (-pkin(1) * t2815 + t2886 * t2922) * t2876;
t2787 = (-pkin(1) * t2814 + t2885 * t2920) * t2874;
t2786 = (-pkin(1) * t2813 + t2884 * t2918) * t2872;
t2776 = t3040 * t3220 + (t2796 * t2921 + t3041 * t3219) * t2904 + t2791 * t3106;
t2775 = t3042 * t3220 + (t2794 * t2919 + t3043 * t3219) * t2904 + t2790 * t3107;
t2774 = t3044 * t3220 + (t2792 * t2917 + t3045 * t3219) * t2904 + t2789 * t3108;
t2773 = -0.2e1 * t2893 * t3041 + 0.2e1 * (-t2791 * t2915 - t2903 * t3040) * t2904 - 0.2e1 * t2791 * t3120;
t2772 = -0.2e1 * t2893 * t3043 + 0.2e1 * (-t2790 * t2913 - t2903 * t3042) * t2904 - 0.2e1 * t2790 * t3121;
t2771 = -0.2e1 * t2893 * t3045 + 0.2e1 * (-t2789 * t2911 - t2903 * t3044) * t2904 - 0.2e1 * t2789 * t3122;
t2760 = t2768 * t3171;
t2759 = t2767 * t3171;
t2758 = t2766 * t3167;
t2757 = t2765 * t3167;
t2756 = t2764 * t3169;
t2755 = t2763 * t3169;
t2751 = t2760 - t3034;
t2750 = -t2760 + 0.2e1 * t3034;
t2749 = t2759 - t3035;
t2748 = -t2759 + 0.2e1 * t3035;
t2747 = t2758 - t3036;
t2746 = -t2758 + 0.2e1 * t3036;
t2745 = t2757 - t3037;
t2744 = -t2757 + 0.2e1 * t3037;
t2743 = t2756 - t3038;
t2742 = -t2756 + 0.2e1 * t3038;
t2741 = t2755 - t3039;
t2740 = -t2755 + 0.2e1 * t3039;
t2739 = t3034 - t2760 / 0.2e1;
t2738 = t3035 - t2759 / 0.2e1;
t2737 = t3036 - t2758 / 0.2e1;
t2736 = t3037 - t2757 / 0.2e1;
t2735 = t3038 - t2756 / 0.2e1;
t2734 = t3039 - t2755 / 0.2e1;
t2733 = (-pkin(1) * t3197 + t2830 * t3145) * t2876;
t2732 = (-pkin(1) * t3199 + t2829 * t3145) * t2876;
t2731 = (-pkin(1) * t3201 + t2828 * t3150) * t2874;
t2730 = (-pkin(1) * t3203 + t2827 * t3150) * t2874;
t2729 = (-pkin(1) * t3205 + t2826 * t3155) * t2872;
t2728 = (-pkin(1) * t3207 + t2825 * t3155) * t2872;
t2724 = t2739 * t3106 + t2750 * t3117 + (t2915 * t2932 + t2938 * t2921) * t2863;
t2723 = t2738 * t3106 + t2748 * t3117 + (t2915 * t2933 + t2939 * t2921) * t2863;
t2722 = t2737 * t3107 + t2746 * t3118 + (t2913 * t2934 + t2940 * t2919) * t2861;
t2721 = t2736 * t3107 + t2744 * t3118 + (t2913 * t2935 + t2941 * t2919) * t2861;
t2720 = t2735 * t3108 + t2742 * t3119 + (t2911 * t2936 + t2942 * t2917) * t2859;
t2719 = t2734 * t3108 + t2740 * t3119 + (t2911 * t2937 + t2943 * t2917) * t2859;
t2718 = -t2739 * t2954 + (-t2938 * t2915 + t2921 * t2932) * t2863;
t2717 = -t2738 * t2954 + (-t2939 * t2915 + t2921 * t2933) * t2863;
t2716 = -t2737 * t2955 + (-t2940 * t2913 + t2919 * t2934) * t2861;
t2715 = -t2736 * t2955 + (-t2941 * t2913 + t2919 * t2935) * t2861;
t2714 = -t2735 * t2956 + (-t2942 * t2911 + t2917 * t2936) * t2859;
t2713 = -t2734 * t2956 + (-t2943 * t2911 + t2917 * t2937) * t2859;
t1 = [(t3085 + t3087 + t3089) * MDP(1) + (qJ(2,1) * t3085 + qJ(2,2) * t3087 + qJ(2,3) * t3089) * t3222 + ((t2733 * t3174 + t2751 * t3197) * t2876 + (t2731 * t3180 + t2747 * t3201) * t2874 + (t2729 * t3186 + t2743 * t3205) * t2872) * MDP(7) + (t2820 * t3072 + t2822 * t3071 + t2824 * t3070) * MDP(8) + (t2810 * t3089 + t2811 * t3087 + t2812 * t3085) * t3221 + (t2843 * t3017 + t2844 * t3012 + t2845 * t3007) * t3110 + (t2840 * t3007 + t2841 * t3017 + t2842 * t3012) * t3109 + (t2860 * t2878 ^ 2 + t2862 * t2879 ^ 2 + t2864 * t2880 ^ 2) * t3208 + ((-t2845 * t2972 + (t2724 * t2876 - t3196) * t2863) * t2830 + (-t2844 * t2976 + (t2722 * t2874 - t3200) * t2861) * t2828 + (-t2843 * t2980 + (t2720 * t2872 - t3204) * t2859) * t2826) * MDP(13) + ((-t2840 * t2972 + (t2718 * t2876 + t3028) * t2863) * t2830 + (-t2842 * t2976 + (t2716 * t2874 + t3030) * t2861) * t2828 + (-t2841 * t2980 + (t2714 * t2872 + t3032) * t2859) * t2826) * MDP(14) + MDP(15) - t3230 * ((-t2750 * t2876 + t3094) * t3174 + (-t2746 * t2874 + t3096) * t3180 + (-t2742 * t2872 + t3098) * t3186); ((t2732 * t3174 + t2749 * t3197) * t2876 + (t2730 * t3180 + t2745 * t3201) * t2874 + (t2728 * t3186 + t2741 * t3205) * t2872) * MDP(7) + ((t2723 * t3173 - t2829 * t3196) * t2863 + (t2721 * t3179 - t2827 * t3200) * t2861 + (t2719 * t3185 - t2825 * t3204) * t2859 + (-t2843 * t2988 - t2844 * t2986 - t2845 * t2984) * t2925) * MDP(13) + ((t2717 * t3173 + t2829 * t3028) * t2863 + (t2715 * t3179 + t2827 * t3030) * t2861 + (t2713 * t3185 + t2825 * t3032) * t2859 + (-t2840 * t2984 - t2841 * t2988 - t2842 * t2986) * t2925) * MDP(14) + t2953 - t3230 * ((-t2740 * t3185 + t2825 * t3098) * t2859 + (-t2744 * t3179 + t2827 * t3096) * t2861 + (-t2748 * t3173 + t2829 * t3094) * t2863); ((t2768 * t3193 + t2788 * t3174) * t2876 + (t2766 * t3194 + t2787 * t3180) * t2874 + (t2764 * t3195 + t2786 * t3186) * t2872) * MDP(7) + (t2774 * t3083 + t2775 * t3081 + t2776 * t3079 - t2764 * t3139 - t2766 * t3136 - t2768 * t3133 + (-t2843 * t2982 - t2844 * t2978 - t2845 * t2974) * t2925) * MDP(13) + (t2764 * t3005 + t2766 * t3004 + t2768 * t3006 + t2771 * t3083 + t2772 * t3081 + t2773 * t3079 + (-t2840 * t2974 - t2841 * t2982 - t2842 * t2978) * t2925) * MDP(14) + t2966 - t3230 * (t2764 * t3075 + t2766 * t3073 + t2768 * t3077 - t2792 * t3083 - t2794 * t3081 - t2796 * t3079); ((t2733 * t3177 + t2751 * t3199) * t2876 + (t2731 * t3183 + t2747 * t3203) * t2874 + (t2729 * t3189 + t2743 * t3207) * t2872) * MDP(7) + ((t2724 * t3176 - t2830 * t3198) * t2863 + (t2722 * t3182 - t2828 * t3202) * t2861 + (t2720 * t3188 - t2826 * t3206) * t2859 + (-t2843 * t2987 - t2844 * t2985 - t2845 * t2983) * t2925) * MDP(13) + ((t2718 * t3176 + t2830 * t3029) * t2863 + (t2716 * t3182 + t2828 * t3031) * t2861 + (t2714 * t3188 + t2826 * t3033) * t2859 + (-t2840 * t2983 - t2841 * t2987 - t2842 * t2985) * t2925) * MDP(14) + t2953 - t3230 * (t2859 * (-t2742 * t3188 + t2826 * t3099) + t2861 * (-t2746 * t3182 + t2828 * t3097) + t2863 * (-t2750 * t3176 + t2830 * t3095)); (t3086 + t3088 + t3090) * MDP(1) + (qJ(2,1) * t3086 + qJ(2,2) * t3088 + qJ(2,3) * t3090) * t3222 + ((t2732 * t3177 + t2749 * t3199) * t2876 + (t2730 * t3183 + t2745 * t3203) * t2874 + (t2728 * t3189 + t2741 * t3207) * t2872) * MDP(7) + (t2819 * t3072 + t2821 * t3071 + t2823 * t3070) * MDP(8) + (t2810 * t3090 + t2811 * t3088 + t2812 * t3086) * t3221 + (t2843 * t3019 + t2844 * t3014 + t2845 * t3009) * t3110 + (t2840 * t3009 + t2841 * t3019 + t2842 * t3014) * t3109 + (t2860 * t2881 ^ 2 + t2862 * t2882 ^ 2 + t2864 * t2883 ^ 2) * t3208 + ((-t2845 * t2971 + (t2723 * t2876 - t3198) * t2863) * t2829 + (-t2844 * t2975 + (t2721 * t2874 - t3202) * t2861) * t2827 + (-t2843 * t2979 + (t2719 * t2872 - t3206) * t2859) * t2825) * MDP(13) + ((-t2840 * t2971 + (t2717 * t2876 + t3029) * t2863) * t2829 + (-t2842 * t2975 + (t2715 * t2874 + t3031) * t2861) * t2827 + (-t2841 * t2979 + (t2713 * t2872 + t3033) * t2859) * t2825) * MDP(14) + MDP(15) - t3230 * ((-t2748 * t2876 + t3095) * t3177 + (-t2744 * t2874 + t3097) * t3183 + (-t2740 * t2872 + t3099) * t3189); ((t2767 * t3193 + t2788 * t3177) * t2876 + (t2765 * t3194 + t2787 * t3183) * t2874 + (t2763 * t3195 + t2786 * t3189) * t2872) * MDP(7) + (t2774 * t3084 + t2775 * t3082 + t2776 * t3080 - t2763 * t3139 - t2765 * t3136 - t2767 * t3133 + (-t2843 * t2981 - t2844 * t2977 - t2845 * t2973) * t2925) * MDP(13) + (t2763 * t3005 + t2765 * t3004 + t2767 * t3006 + t2771 * t3084 + t2772 * t3082 + t2773 * t3080 + (-t2840 * t2973 - t2841 * t2981 - t2842 * t2977) * t2925) * MDP(14) + t2967 - t3230 * (t2763 * t3075 + t2765 * t3073 + t2767 * t3077 - t2792 * t3084 - t2794 * t3082 - t2796 * t3080); ((t2733 * t2922 + t2751 * t2815) * t2876 + (t2731 * t2920 + t2747 * t2814) * t2874 + (t2729 * t2918 + t2743 * t2813) * t2872) * MDP(7) + (t2720 * t3141 + t2722 * t3138 + t2724 * t3135 - t2826 * t3027 - t2828 * t3025 - t2830 * t3023) * MDP(13) + (t2714 * t3141 + t2716 * t3138 + t2718 * t3135 + t2826 * t3026 + t2828 * t3024 + t2830 * t3022) * MDP(14) + t2966 - t3230 * (-t2742 * t3141 - t2746 * t3138 - t2750 * t3135 + t2826 * t3093 + t2828 * t3092 + t2830 * t3091); ((t2732 * t2922 + t2749 * t2815) * t2876 + (t2730 * t2920 + t2745 * t2814) * t2874 + (t2728 * t2918 + t2741 * t2813) * t2872) * MDP(7) + (t2719 * t3141 + t2721 * t3138 + t2723 * t3135 - t2825 * t3027 - t2827 * t3025 - t2829 * t3023) * MDP(13) + (t2713 * t3141 + t2715 * t3138 + t2717 * t3135 + t2825 * t3026 + t2827 * t3024 + t2829 * t3022) * MDP(14) + t2967 - t3230 * (-t2740 * t3141 - t2744 * t3138 - t2748 * t3135 + t2825 * t3093 + t2827 * t3092 + t2829 * t3091); (t3134 + t3137 + t3140) * MDP(1) + (qJ(2,1) * t3134 + qJ(2,2) * t3137 + qJ(2,3) * t3140) * t3222 + (t2793 * t2804 + t2795 * t2805 + t2797 * t2806) * MDP(7) + (t2837 * t3140 + t2838 * t3137 + t2839 * t3134) * MDP(8) + (t2810 * t3140 + t2811 * t3137 + t2812 * t3134) * t3221 + MDP(15) + (t2788 * t2876 * MDP(7) + (t2776 * t2876 - t2840 * t3190) * MDP(13) + (t2773 * t2876 + t2815 * t3163) * MDP(14)) * t2922 + (t2787 * t2874 * MDP(7) + (t2775 * t2874 - t2842 * t3191) * MDP(13) + (t2772 * t2874 + t2814 * t3164) * MDP(14)) * t2920 + (t2786 * t2872 * MDP(7) + (t2774 * t2872 - t2841 * t3192) * MDP(13) + (t2771 * t2872 + t2813 * t3165) * MDP(14)) * t2918 + t3230 * (t2918 * (t2792 * t2872 - t3192) + t2920 * (t2794 * t2874 - t3191) + t2922 * (t2796 * t2876 - t3190));];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
