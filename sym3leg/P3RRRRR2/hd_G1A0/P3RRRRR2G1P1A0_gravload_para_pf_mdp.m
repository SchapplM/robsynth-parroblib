% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G1P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR2G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:05:09
% EndTime: 2020-03-09 21:05:10
% DurationCPUTime: 1.20s
% Computational Cost: add. (918->159), mult. (876->250), div. (255->14), fcn. (927->60), ass. (0->133)
t2843 = legFrame(3,3);
t2825 = sin(t2843);
t2904 = cos(t2843);
t2786 = t2825 * g(1) - t2904 * g(2);
t2789 = t2904 * g(1) + t2825 * g(2);
t2840 = qJ(1,3) + qJ(2,3);
t2819 = sin(t2840);
t2822 = cos(t2840);
t2768 = t2786 * t2822 + t2789 * t2819;
t2846 = sin(qJ(3,3));
t2908 = t2768 * t2846 ^ 2;
t2844 = legFrame(2,3);
t2826 = sin(t2844);
t2903 = cos(t2844);
t2787 = t2826 * g(1) - t2903 * g(2);
t2790 = t2903 * g(1) + t2826 * g(2);
t2841 = qJ(1,2) + qJ(2,2);
t2820 = sin(t2841);
t2823 = cos(t2841);
t2769 = t2787 * t2823 + t2790 * t2820;
t2848 = sin(qJ(3,2));
t2907 = t2769 * t2848 ^ 2;
t2845 = legFrame(1,3);
t2827 = sin(t2845);
t2902 = cos(t2845);
t2788 = t2827 * g(1) - t2902 * g(2);
t2791 = t2902 * g(1) + t2827 * g(2);
t2842 = qJ(1,1) + qJ(2,1);
t2821 = sin(t2842);
t2824 = cos(t2842);
t2770 = t2788 * t2824 + t2791 * t2821;
t2850 = sin(qJ(3,1));
t2906 = t2770 * t2850 ^ 2;
t2905 = -2 * pkin(1);
t2818 = qJ(1,1) + t2845;
t2817 = qJ(1,2) + t2844;
t2816 = qJ(1,3) + t2843;
t2852 = cos(qJ(3,3));
t2901 = t2768 * t2852;
t2854 = cos(qJ(3,2));
t2900 = t2769 * t2854;
t2856 = cos(qJ(3,1));
t2899 = t2770 * t2856;
t2810 = qJ(2,3) + t2816;
t2804 = qJ(3,3) + t2810;
t2805 = -qJ(3,3) + t2810;
t2780 = sin(t2816) * t2905 + (-sin(t2805) - sin(t2804)) * pkin(2);
t2792 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t2898 = t2780 * t2792;
t2811 = qJ(2,2) + t2817;
t2806 = qJ(3,2) + t2811;
t2807 = -qJ(3,2) + t2811;
t2781 = sin(t2817) * t2905 + (-sin(t2807) - sin(t2806)) * pkin(2);
t2793 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t2897 = t2781 * t2793;
t2812 = qJ(2,1) + t2818;
t2808 = qJ(3,1) + t2812;
t2809 = -qJ(3,1) + t2812;
t2782 = sin(t2818) * t2905 + (-sin(t2809) - sin(t2808)) * pkin(2);
t2794 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t2896 = t2782 * t2794;
t2783 = cos(t2816) * t2905 + (-cos(t2805) - cos(t2804)) * pkin(2);
t2895 = t2783 * t2792;
t2784 = cos(t2817) * t2905 + (-cos(t2807) - cos(t2806)) * pkin(2);
t2894 = t2784 * t2793;
t2785 = cos(t2818) * t2905 + (-cos(t2809) - cos(t2808)) * pkin(2);
t2893 = t2785 * t2794;
t2798 = sin(t2810);
t2829 = 0.1e1 / sin(qJ(2,3));
t2892 = t2798 * t2829;
t2799 = sin(t2811);
t2831 = 0.1e1 / sin(qJ(2,2));
t2891 = t2799 * t2831;
t2800 = sin(t2812);
t2833 = 0.1e1 / sin(qJ(2,1));
t2890 = t2800 * t2833;
t2801 = cos(t2810);
t2889 = t2801 * t2829;
t2802 = cos(t2811);
t2888 = t2802 * t2831;
t2803 = cos(t2812);
t2887 = t2803 * t2833;
t2886 = t2829 * t2846;
t2885 = t2831 * t2848;
t2884 = t2833 * t2850;
t2883 = t2768 * t2792 * t2846;
t2882 = t2792 * t2901;
t2881 = t2768 * t2886;
t2880 = t2829 * t2901;
t2879 = t2769 * t2793 * t2848;
t2878 = t2793 * t2900;
t2877 = t2769 * t2885;
t2876 = t2831 * t2900;
t2875 = t2770 * t2794 * t2850;
t2874 = t2794 * t2899;
t2873 = t2770 * t2884;
t2872 = t2833 * t2899;
t2795 = cos(qJ(2,3)) * pkin(1) + t2852 * pkin(2);
t2871 = t2795 * t2829 / t2852 ^ 2;
t2796 = cos(qJ(2,2)) * pkin(1) + t2854 * pkin(2);
t2870 = t2796 * t2831 / t2854 ^ 2;
t2797 = cos(qJ(2,1)) * pkin(1) + t2856 * pkin(2);
t2869 = t2797 * t2833 / t2856 ^ 2;
t2834 = 0.1e1 / t2852;
t2868 = t2834 * t2886;
t2836 = 0.1e1 / t2854;
t2867 = t2836 * t2885;
t2838 = 0.1e1 / t2856;
t2866 = t2838 * t2884;
t2865 = t2768 * t2868;
t2864 = t2769 * t2867;
t2863 = t2770 * t2866;
t2862 = t2846 * t2871;
t2861 = t2848 * t2870;
t2860 = t2850 * t2869;
t2771 = -t2786 * t2819 + t2789 * t2822;
t2772 = -t2787 * t2820 + t2790 * t2823;
t2773 = -t2788 * t2821 + t2791 * t2824;
t2859 = 1 / pkin(1);
t2858 = 0.1e1 / pkin(2);
t2857 = cos(qJ(1,1));
t2855 = cos(qJ(1,2));
t2853 = cos(qJ(1,3));
t2851 = sin(qJ(1,1));
t2849 = sin(qJ(1,2));
t2847 = sin(qJ(1,3));
t2779 = -t2788 * t2851 + t2791 * t2857;
t2778 = -t2787 * t2849 + t2790 * t2855;
t2777 = -t2786 * t2847 + t2789 * t2853;
t2776 = t2788 * t2857 + t2791 * t2851;
t2775 = t2787 * t2855 + t2790 * t2849;
t2774 = t2786 * t2853 + t2789 * t2847;
t1 = [-g(1) * MDP(14) + ((t2774 * t2889 + t2775 * t2888 + t2776 * t2887) * MDP(2) + (t2777 * t2889 + t2778 * t2888 + t2779 * t2887) * MDP(3) + (t2768 * t2889 + t2769 * t2888 + t2770 * t2887) * MDP(5) + (t2771 * t2889 + t2772 * t2888 + t2773 * t2887) * MDP(6) + (t2801 * t2880 + t2802 * t2876 + t2803 * t2872) * MDP(12) + (-t2801 * t2881 - t2802 * t2877 - t2803 * t2873) * MDP(13) + ((t2768 * t2895 + t2769 * t2894 + t2770 * t2893) * MDP(5) + (t2771 * t2895 + t2772 * t2894 + t2773 * t2893) * MDP(6) + (t2783 * t2882 + t2784 * t2878 + t2785 * t2874) * MDP(12) + (-t2783 * t2883 - t2784 * t2879 - t2785 * t2875) * MDP(13)) * t2858) * t2859; -g(2) * MDP(14) + ((t2774 * t2892 + t2775 * t2891 + t2776 * t2890) * MDP(2) + (t2777 * t2892 + t2778 * t2891 + t2779 * t2890) * MDP(3) + (t2768 * t2892 + t2769 * t2891 + t2770 * t2890) * MDP(5) + (t2771 * t2892 + t2772 * t2891 + t2773 * t2890) * MDP(6) + (t2798 * t2880 + t2799 * t2876 + t2800 * t2872) * MDP(12) + (-t2798 * t2881 - t2799 * t2877 - t2800 * t2873) * MDP(13) + ((t2768 * t2898 + t2769 * t2897 + t2770 * t2896) * MDP(5) + (t2771 * t2898 + t2772 * t2897 + t2773 * t2896) * MDP(6) + (t2780 * t2882 + t2781 * t2878 + t2782 * t2874) * MDP(12) + (-t2780 * t2883 - t2781 * t2879 - t2782 * t2875) * MDP(13)) * t2858) * t2859; -g(3) * MDP(14) + ((t2838 * (-g(3) * t2856 + t2773 * t2850) + t2836 * (-g(3) * t2854 + t2772 * t2848) + t2834 * (-g(3) * t2852 + t2771 * t2846)) * MDP(12) + (t2838 * (g(3) * t2850 + t2773 * t2856) + t2836 * (g(3) * t2848 + t2772 * t2854) + t2834 * (g(3) * t2846 + t2771 * t2852)) * MDP(13)) * t2858 + ((t2774 * t2868 + t2775 * t2867 + t2776 * t2866) * MDP(2) + (t2777 * t2868 + t2778 * t2867 + t2779 * t2866) * MDP(3) + (t2863 + t2864 + t2865) * MDP(5) + (t2771 * t2868 + t2772 * t2867 + t2773 * t2866) * MDP(6) + (t2873 + t2877 + t2881) * MDP(12) + (-t2834 * t2829 * t2908 - t2836 * t2831 * t2907 - t2838 * t2833 * t2906) * MDP(13) + ((-t2768 * t2862 - t2769 * t2861 - t2770 * t2860) * MDP(5) + (-t2771 * t2862 - t2772 * t2861 - t2773 * t2860) * MDP(6) + (-t2795 * t2865 - t2796 * t2864 - t2797 * t2863) * MDP(12) + (t2869 * t2906 + t2870 * t2907 + t2871 * t2908) * MDP(13)) * t2858) * t2859;];
taugX  = t1;
