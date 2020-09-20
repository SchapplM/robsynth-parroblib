% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR2G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G3P1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G3P1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:09
% EndTime: 2020-03-09 21:20:10
% DurationCPUTime: 0.90s
% Computational Cost: add. (4402->127), mult. (3566->250), div. (1356->9), fcn. (2406->18), ass. (0->132)
t901 = sin(qJ(3,3));
t890 = 0.1e1 / t901 ^ 2;
t904 = cos(qJ(3,3));
t974 = t890 * t904;
t902 = sin(qJ(3,2));
t893 = 0.1e1 / t902 ^ 2;
t905 = cos(qJ(3,2));
t973 = t893 * t905;
t903 = sin(qJ(3,1));
t896 = 0.1e1 / t903 ^ 2;
t906 = cos(qJ(3,1));
t972 = t896 * t906;
t971 = 2 * pkin(2);
t889 = 0.1e1 / t901;
t892 = 0.1e1 / t902;
t895 = 0.1e1 / t903;
t912 = 0.1e1 / pkin(1) ^ 2;
t945 = t890 * t912;
t910 = 1 / pkin(2);
t911 = 1 / pkin(1);
t940 = t910 * t911;
t886 = -legFrame(3,2) + qJ(2,3);
t874 = sin(t886);
t877 = cos(t886);
t907 = xDP(2);
t908 = xDP(1);
t880 = qJ(3,3) + t886;
t868 = sin(t880);
t871 = cos(t880);
t921 = t868 * t908 - t871 * t907;
t850 = t921 * pkin(2) + (t874 * t908 - t877 * t907) * pkin(1);
t961 = t850 * t889;
t847 = t940 * t961;
t955 = t921 * t889;
t853 = t911 * t955;
t842 = -t853 + t847;
t964 = t842 * t850;
t835 = t945 * t964;
t946 = t889 * t904;
t838 = -pkin(2) * t842 + t921 * t946;
t926 = (pkin(1) * t904 + pkin(2)) * t910 * t964;
t909 = pkin(2) ^ 2;
t967 = (t842 * t909 + ((-t853 + t847 / 0.2e1) * t904 * t971 - t955) * pkin(1)) * t910;
t826 = t835 + (-t926 - (-t838 - t967) * t921) * t945;
t970 = t826 * t889;
t943 = t893 * t912;
t887 = -legFrame(2,2) + qJ(2,2);
t875 = sin(t887);
t878 = cos(t887);
t881 = qJ(3,2) + t887;
t869 = sin(t881);
t872 = cos(t881);
t920 = t869 * t908 - t872 * t907;
t851 = t920 * pkin(2) + (t875 * t908 - t878 * t907) * pkin(1);
t960 = t851 * t892;
t848 = t940 * t960;
t954 = t920 * t892;
t854 = t911 * t954;
t844 = -t854 + t848;
t963 = t844 * t851;
t836 = t943 * t963;
t944 = t892 * t905;
t839 = -pkin(2) * t844 + t920 * t944;
t924 = (pkin(1) * t905 + pkin(2)) * t910 * t963;
t966 = (t844 * t909 + ((-t854 + t848 / 0.2e1) * t905 * t971 - t954) * pkin(1)) * t910;
t827 = t836 + (-t924 - (-t839 - t966) * t920) * t943;
t969 = t827 * t892;
t941 = t896 * t912;
t888 = -legFrame(1,2) + qJ(2,1);
t876 = sin(t888);
t879 = cos(t888);
t882 = qJ(3,1) + t888;
t870 = sin(t882);
t873 = cos(t882);
t919 = t870 * t908 - t873 * t907;
t852 = t919 * pkin(2) + (t876 * t908 - t879 * t907) * pkin(1);
t959 = t852 * t895;
t849 = t940 * t959;
t953 = t919 * t895;
t855 = t911 * t953;
t846 = -t855 + t849;
t962 = t846 * t852;
t837 = t941 * t962;
t942 = t895 * t906;
t840 = -pkin(2) * t846 + t919 * t942;
t922 = (pkin(1) * t906 + pkin(2)) * t910 * t962;
t965 = (t846 * t909 + ((-t855 + t849 / 0.2e1) * t906 * t971 - t953) * pkin(1)) * t910;
t828 = t837 + (-t922 - (-t840 - t965) * t919) * t941;
t968 = t828 * t895;
t856 = t921 ^ 2;
t958 = t856 * t890;
t857 = t920 ^ 2;
t957 = t857 * t893;
t858 = t919 ^ 2;
t956 = t858 * t896;
t952 = t868 * t889;
t951 = t869 * t892;
t950 = t870 * t895;
t949 = t871 * t889;
t948 = t872 * t892;
t947 = t873 * t895;
t829 = t838 * t921 * t945 + t835;
t939 = t829 * t946;
t830 = t839 * t920 * t943 + t836;
t938 = t830 * t944;
t831 = t840 * t919 * t941 + t837;
t937 = t831 * t942;
t841 = -0.2e1 * t853 + t847;
t936 = t841 * t961;
t843 = -0.2e1 * t854 + t848;
t935 = t843 * t960;
t845 = -0.2e1 * t855 + t849;
t934 = t845 * t959;
t933 = t856 * t889 * t974;
t932 = t857 * t892 * t973;
t931 = t858 * t895 * t972;
t823 = 0.2e1 * t835 + (-t926 - (-0.2e1 * t838 - t967) * t921) * t945;
t930 = t823 * t946;
t824 = 0.2e1 * t836 + (-t924 - (-0.2e1 * t839 - t966) * t920) * t943;
t929 = t824 * t944;
t825 = 0.2e1 * t837 + (-t922 - (-0.2e1 * t840 - t965) * t919) * t941;
t928 = t825 * t942;
t927 = t841 * t850 * t974;
t925 = t843 * t851 * t973;
t923 = t845 * t852 * t972;
t867 = -pkin(1) * t879 - pkin(2) * t873;
t866 = -pkin(1) * t878 - pkin(2) * t872;
t865 = -pkin(1) * t877 - pkin(2) * t871;
t864 = pkin(1) * t876 + pkin(2) * t870;
t863 = pkin(1) * t875 + pkin(2) * t869;
t862 = pkin(1) * t874 + pkin(2) * t868;
t1 = [(-t868 * t930 - t869 * t929 - t870 * t928) * MDP(6) + (t868 * t823 + t869 * t824 + t870 * t825) * MDP(7) + ((-t829 * t952 - t830 * t951 - t831 * t950) * MDP(2) + (-t826 * t952 - t827 * t951 - t828 * t950) * MDP(5)) * t911 + ((t862 * t939 + t863 * t938 + t864 * t937) * MDP(6) + (-t862 * t829 - t863 * t830 - t864 * t831) * MDP(7) + ((t862 * t958 + t863 * t957 + t864 * t956) * MDP(6) + (t862 * t933 + t863 * t932 + t864 * t931) * MDP(7)) * t912 + ((t862 * t970 + t863 * t969 + t864 * t968) * MDP(5) + (t868 * t936 + t869 * t935 + t870 * t934) * MDP(6) + (t868 * t927 + t869 * t925 + t870 * t923) * MDP(7)) * t911) * t910; (t871 * t930 + t872 * t929 + t873 * t928) * MDP(6) + (-t871 * t823 - t872 * t824 - t873 * t825) * MDP(7) + ((t829 * t949 + t830 * t948 + t831 * t947) * MDP(2) + (t826 * t949 + t827 * t948 + t828 * t947) * MDP(5)) * t911 + ((t865 * t939 + t866 * t938 + t867 * t937) * MDP(6) + (-t865 * t829 - t866 * t830 - t867 * t831) * MDP(7) + ((t865 * t958 + t866 * t957 + t867 * t956) * MDP(6) + (t865 * t933 + t866 * t932 + t867 * t931) * MDP(7)) * t912 + ((t865 * t970 + t866 * t969 + t867 * t968) * MDP(5) + (-t871 * t936 - t872 * t935 - t873 * t934) * MDP(6) + (-t871 * t927 - t872 * t925 - t873 * t923) * MDP(7)) * t911) * t910; 0;];
taucX  = t1;
