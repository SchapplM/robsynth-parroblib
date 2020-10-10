% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:29
% EndTime: 2020-03-09 21:02:31
% DurationCPUTime: 2.17s
% Computational Cost: add. (1965->261), mult. (5692->512), div. (2718->17), fcn. (5772->18), ass. (0->193)
t866 = legFrame(3,2);
t838 = sin(t866);
t841 = cos(t866);
t872 = sin(qJ(3,3));
t878 = cos(qJ(3,3));
t885 = xDP(2);
t886 = xDP(1);
t879 = cos(qJ(2,3));
t884 = xDP(3);
t949 = t879 * t884;
t820 = -t872 * t949 + (t838 * t885 - t841 * t886) * t878;
t817 = t820 ^ 2;
t873 = sin(qJ(2,3));
t846 = 0.1e1 / t873 ^ 2;
t983 = t817 * t846;
t867 = legFrame(2,2);
t839 = sin(t867);
t842 = cos(t867);
t874 = sin(qJ(3,2));
t880 = cos(qJ(3,2));
t881 = cos(qJ(2,2));
t947 = t881 * t884;
t821 = -t874 * t947 + (t839 * t885 - t842 * t886) * t880;
t818 = t821 ^ 2;
t875 = sin(qJ(2,2));
t849 = 0.1e1 / t875 ^ 2;
t982 = t818 * t849;
t868 = legFrame(1,2);
t840 = sin(t868);
t843 = cos(t868);
t876 = sin(qJ(3,1));
t882 = cos(qJ(3,1));
t883 = cos(qJ(2,1));
t945 = t883 * t884;
t822 = -t876 * t945 + (t840 * t885 - t843 * t886) * t882;
t819 = t822 ^ 2;
t877 = sin(qJ(2,1));
t852 = 0.1e1 / t877 ^ 2;
t981 = t819 * t852;
t887 = 0.1e1 / pkin(2);
t845 = 0.1e1 / t873;
t848 = 0.1e1 / t875;
t851 = 0.1e1 / t877;
t853 = 0.1e1 / t878;
t857 = 0.1e1 / t880;
t861 = 0.1e1 / t882;
t888 = 0.1e1 / pkin(2) ^ 2;
t898 = t878 ^ 2;
t854 = 0.1e1 / t898;
t901 = t880 ^ 2;
t858 = 0.1e1 / t901;
t904 = t882 ^ 2;
t862 = 0.1e1 / t904;
t987 = 2 * MDP(6);
t829 = -t838 * t879 + t841 * t873;
t830 = t873 * t838 + t841 * t879;
t855 = t853 * t854;
t869 = xDDP(3);
t870 = xDDP(2);
t871 = xDDP(1);
t865 = t884 ^ 2;
t953 = t865 * t887;
t962 = t853 * t872;
t808 = -t838 * g(1) - t841 * g(2) + (t829 * t870 + t830 * t871 + t869 * t962 + (t887 * t983 + t953) * t855) * t845;
t835 = t841 * g(1) - t838 * g(2);
t802 = t808 * t879 + t835 * t873;
t986 = t802 * t845;
t831 = -t839 * t881 + t842 * t875;
t832 = t875 * t839 + t842 * t881;
t859 = t857 * t858;
t959 = t857 * t874;
t809 = -t839 * g(1) - t842 * g(2) + (t831 * t870 + t832 * t871 + t869 * t959 + (t887 * t982 + t953) * t859) * t848;
t836 = t842 * g(1) - t839 * g(2);
t803 = t809 * t881 + t836 * t875;
t985 = t803 * t848;
t833 = -t840 * t883 + t843 * t877;
t834 = t877 * t840 + t843 * t883;
t863 = t861 * t862;
t956 = t861 * t876;
t810 = -t840 * g(1) - t843 * g(2) + (t833 * t870 + t834 * t871 + t869 * t956 + (t887 * t981 + t953) * t863) * t851;
t837 = t843 * g(1) - t840 * g(2);
t804 = t810 * t883 + t837 * t877;
t984 = t804 * t851;
t952 = t865 * t888;
t921 = t872 * t952;
t951 = t869 * t887;
t826 = t853 * t951 + t855 * t921;
t980 = t826 * t872;
t979 = t826 * t878;
t920 = t874 * t952;
t827 = t857 * t951 + t859 * t920;
t978 = t827 * t874;
t977 = t827 * t880;
t919 = t876 * t952;
t828 = t861 * t951 + t863 * t919;
t976 = t828 * t876;
t975 = t828 * t882;
t974 = t829 * t845;
t973 = t830 * t845;
t972 = t831 * t848;
t971 = t832 * t848;
t970 = t833 * t851;
t969 = t834 * t851;
t968 = t845 * t853;
t967 = t845 * t879;
t966 = t848 * t857;
t965 = t848 * t881;
t964 = t851 * t861;
t963 = t851 * t883;
t961 = t854 * t879;
t856 = 0.1e1 / t898 ^ 2;
t960 = t856 * t888;
t958 = t858 * t881;
t860 = 0.1e1 / t901 ^ 2;
t957 = t860 * t888;
t955 = t862 * t883;
t864 = 0.1e1 / t904 ^ 2;
t954 = t864 * t888;
t924 = t872 * t961;
t931 = t846 * t960;
t944 = t884 * t887;
t799 = -(-t872 * t873 * t884 + t820 * t967) * t820 * t931 + (-t869 * t924 + (t838 * t870 - t841 * t871 - (-t820 * t872 + t949) * t855 * t944) * t853) * t845 * t887;
t950 = t878 * t799;
t923 = t874 * t958;
t928 = t849 * t957;
t800 = -(-t874 * t875 * t884 + t821 * t965) * t821 * t928 + (-t869 * t923 + (t839 * t870 - t842 * t871 - (-t821 * t874 + t947) * t859 * t944) * t857) * t848 * t887;
t948 = t880 * t800;
t922 = t876 * t955;
t925 = t852 * t954;
t801 = -(-t876 * t877 * t884 + t822 * t963) * t822 * t925 + (-t869 * t922 + (t840 * t870 - t843 * t871 - (-t822 * t876 + t945) * t863 * t944) * t861) * t851 * t887;
t946 = t882 * t801;
t943 = t884 * t888;
t942 = t856 * t983;
t941 = t860 * t982;
t940 = t864 * t981;
t939 = t838 * t968;
t938 = t839 * t966;
t937 = t840 * t964;
t936 = t841 * t968;
t935 = t842 * t966;
t934 = t843 * t964;
t933 = t845 * t962;
t932 = t845 * t961;
t930 = t848 * t959;
t929 = t848 * t958;
t927 = t851 * t956;
t926 = t851 * t955;
t918 = t802 * t933;
t917 = t803 * t930;
t916 = t804 * t927;
t915 = t820 * t845 * t943;
t914 = t821 * t848 * t943;
t913 = t822 * t851 * t943;
t912 = t845 * t924;
t911 = t848 * t923;
t910 = t851 * t922;
t909 = 0.2e1 * t854 * t915;
t908 = 0.2e1 * t858 * t914;
t907 = 0.2e1 * t862 * t913;
t850 = t876 ^ 2;
t847 = t874 ^ 2;
t844 = t872 ^ 2;
t825 = t861 * t952 + t976;
t824 = t857 * t952 + t978;
t823 = t853 * t952 + t980;
t816 = -t862 * t919 + t975;
t815 = -t858 * t920 + t977;
t814 = -t854 * t921 + t979;
t813 = (t862 * t865 + t940) * t888;
t812 = (t858 * t865 + t941) * t888;
t811 = (t854 * t865 + t942) * t888;
t807 = -t810 * t877 + t837 * t883;
t806 = -t809 * t875 + t836 * t881;
t805 = -t808 * t873 + t835 * t879;
t798 = -t819 * t851 * t954 + t883 * t801;
t797 = -t818 * t848 * t957 + t881 * t800;
t796 = -t817 * t845 * t960 + t879 * t799;
t795 = -t819 * t883 * t925 - t801 * t877;
t794 = -t818 * t881 * t928 - t800 * t875;
t793 = -t817 * t879 * t931 - t799 * t873;
t792 = t850 * t801 + t876 * t907;
t791 = t847 * t800 + t874 * t908;
t790 = t844 * t799 + t872 * t909;
t789 = t876 * t946 + (0.2e1 * t861 - t863) * t913;
t788 = t874 * t948 + (0.2e1 * t857 - t859) * t914;
t787 = t872 * t950 + (0.2e1 * t853 - t855) * t915;
t786 = (t813 * t876 - t975) * t877 - t883 * (t876 * t801 + t907);
t785 = (t812 * t874 - t977) * t875 - t881 * (t874 * t800 + t908);
t784 = (t811 * t872 - t979) * t873 - t879 * (t872 * t799 + t909);
t783 = (-t813 * t882 - t976) * t877 + (-0.2e1 * t876 * t863 * t913 + t946) * t883;
t782 = (-t812 * t880 - t978) * t875 + (-0.2e1 * t874 * t859 * t914 + t948) * t881;
t781 = (-t811 * t878 - t980) * t873 + (-0.2e1 * t872 * t855 * t915 + t950) * t879;
t1 = [(t808 * t973 + t809 * t971 + t810 * t969) * MDP(1) + (t796 * t973 + t797 * t971 + t798 * t969) * MDP(3) + (t793 * t973 + t794 * t971 + t795 * t969) * MDP(4) + (t781 * t973 + t782 * t971 + t783 * t969) * MDP(10) + (t784 * t973 + t785 * t971 + t786 * t969) * MDP(11) + (t871 - g(1)) * MDP(12) + ((-t799 * t936 - t800 * t935 - t801 * t934) * MDP(2) + (-t802 * t936 - t803 * t935 - t804 * t934) * MDP(3) + (-t805 * t936 - t806 * t935 - t807 * t934) * MDP(4) + (-t790 * t936 - t791 * t935 - t792 * t934) * MDP(5) + (-t823 * t936 - t824 * t935 - t825 * t934) * MDP(7) + (-t814 * t936 - t815 * t935 - t816 * t934) * MDP(8) + (-t841 * t986 - t842 * t985 - t843 * t984) * MDP(10) + (t841 * t918 + t842 * t917 + t843 * t916) * MDP(11) + (-t787 * t936 - t788 * t935 - t789 * t934) * t987) * t887; (t808 * t974 + t809 * t972 + t810 * t970) * MDP(1) + (t796 * t974 + t797 * t972 + t798 * t970) * MDP(3) + (t793 * t974 + t794 * t972 + t795 * t970) * MDP(4) + (t781 * t974 + t782 * t972 + t783 * t970) * MDP(10) + (t784 * t974 + t785 * t972 + t786 * t970) * MDP(11) + (t870 - g(2)) * MDP(12) + ((t799 * t939 + t800 * t938 + t801 * t937) * MDP(2) + (t802 * t939 + t803 * t938 + t804 * t937) * MDP(3) + (t805 * t939 + t806 * t938 + t807 * t937) * MDP(4) + (t790 * t939 + t791 * t938 + t792 * t937) * MDP(5) + (t823 * t939 + t824 * t938 + t825 * t937) * MDP(7) + (t814 * t939 + t815 * t938 + t816 * t937) * MDP(8) + (t838 * t986 + t839 * t985 + t840 * t984) * MDP(10) + (-t838 * t918 - t839 * t917 - t840 * t916) * MDP(11) + (t787 * t939 + t788 * t938 + t789 * t937) * t987) * t887; (t808 * t933 + t809 * t930 + t810 * t927) * MDP(1) + (t796 * t933 + t797 * t930 + t798 * t927) * MDP(3) + (t793 * t933 + t794 * t930 + t795 * t927) * MDP(4) + (t781 * t933 + t782 * t930 + t783 * t927) * MDP(10) + (t784 * t933 + t785 * t930 + t786 * t927) * MDP(11) + (t869 - g(3)) * MDP(12) + ((-t799 * t912 - t800 * t911 - t801 * t910) * MDP(2) + (-t802 * t912 - t803 * t911 - t804 * t910) * MDP(3) + (-t805 * t912 - t806 * t911 - t807 * t910) * MDP(4) + ((-t872 * t942 - t874 * t941 - t876 * t940) * t888 - t790 * t912 - t791 * t911 - t792 * t910) * MDP(5) + (-0.2e1 * t787 * t912 - 0.2e1 * t788 * t911 - 0.2e1 * t789 * t910 + (t861 * (-0.2e1 * t862 + t864) * t981 + t857 * (-0.2e1 * t858 + t860) * t982 + t853 * (-0.2e1 * t854 + t856) * t983) * t888) * MDP(6) + ((t801 * t861 - t825 * t926) * t876 + (t800 * t857 - t824 * t929) * t874 + (t799 * t853 - t823 * t932) * t872) * MDP(7) + (-t814 * t912 - t815 * t911 - t816 * t910 + t799 + t800 + t801) * MDP(8) + (t826 * t853 + t827 * t857 + t828 * t861) * MDP(9) + ((-g(3) * t882 + (-t804 * t963 + t807) * t876) * t861 + (-g(3) * t880 + (-t803 * t965 + t806) * t874) * t857 + (-g(3) * t878 + (-t802 * t967 + t805) * t872) * t853) * MDP(10) + (t850 * t804 * t926 + t861 * (g(3) * t876 + t807 * t882) + t847 * t803 * t929 + t857 * (g(3) * t874 + t806 * t880) + t844 * t802 * t932 + t853 * (g(3) * t872 + t805 * t878)) * MDP(11)) * t887;];
tauX  = t1;
