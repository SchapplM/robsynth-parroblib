% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x10]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRR1G1P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:37
% EndTime: 2019-05-03 15:38:38
% DurationCPUTime: 0.52s
% Computational Cost: add. (788->121), mult. (988->239), div. (126->5), fcn. (1084->20), ass. (0->104)
t874 = qJ(1,3) + qJ(2,3);
t860 = sin(t874);
t863 = cos(t874);
t880 = sin(qJ(1,3));
t883 = cos(qJ(1,3));
t912 = 0.1e1 / (t860 * t883 - t863 * t880);
t875 = qJ(1,2) + qJ(2,2);
t861 = sin(t875);
t864 = cos(t875);
t881 = sin(qJ(1,2));
t884 = cos(qJ(1,2));
t911 = 0.1e1 / (t861 * t884 - t864 * t881);
t876 = qJ(1,1) + qJ(2,1);
t862 = sin(t876);
t865 = cos(t876);
t882 = sin(qJ(1,1));
t885 = cos(qJ(1,1));
t910 = 0.1e1 / (t862 * t885 - t865 * t882);
t886 = xP(3);
t872 = sin(t886);
t873 = cos(t886);
t887 = koppelP(3,2);
t890 = koppelP(3,1);
t837 = t872 * t890 + t873 * t887;
t840 = -t872 * t887 + t873 * t890;
t877 = legFrame(3,3);
t866 = sin(t877);
t869 = cos(t877);
t795 = (t837 * t869 - t840 * t866) * t863 - (t837 * t866 + t840 * t869) * t860;
t909 = t795 * t912;
t888 = koppelP(2,2);
t891 = koppelP(2,1);
t838 = t872 * t891 + t873 * t888;
t841 = -t872 * t888 + t873 * t891;
t878 = legFrame(2,3);
t867 = sin(t878);
t870 = cos(t878);
t796 = (t838 * t870 - t841 * t867) * t864 - (t838 * t867 + t841 * t870) * t861;
t908 = t796 * t911;
t889 = koppelP(1,2);
t892 = koppelP(1,1);
t839 = t872 * t892 + t873 * t889;
t842 = -t872 * t889 + t873 * t892;
t879 = legFrame(1,3);
t868 = sin(t879);
t871 = cos(t879);
t797 = (t839 * t871 - t842 * t868) * t865 - (t839 * t868 + t842 * t871) * t862;
t907 = t797 * t910;
t843 = -t866 * g(1) + g(2) * t869;
t846 = g(1) * t869 + g(2) * t866;
t810 = -t843 * t863 + t846 * t860;
t906 = t810 * t912;
t811 = t843 * t860 + t846 * t863;
t905 = t811 * t912;
t844 = -t867 * g(1) + g(2) * t870;
t847 = g(1) * t870 + g(2) * t867;
t812 = -t844 * t864 + t847 * t861;
t904 = t812 * t911;
t813 = t844 * t861 + t847 * t864;
t903 = t813 * t911;
t845 = -t868 * g(1) + g(2) * t871;
t848 = g(1) * t871 + g(2) * t868;
t814 = -t845 * t865 + t848 * t862;
t902 = t814 * t910;
t815 = t845 * t862 + t848 * t865;
t901 = t815 * t910;
t822 = t860 * t869 + t863 * t866;
t900 = t822 * t912;
t823 = t860 * t866 - t863 * t869;
t899 = t823 * t912;
t824 = t861 * t870 + t864 * t867;
t898 = t824 * t911;
t825 = t861 * t867 - t864 * t870;
t897 = t825 * t911;
t826 = t862 * t871 + t865 * t868;
t896 = t826 * t910;
t827 = t862 * t868 - t865 * t871;
t895 = t827 * t910;
t894 = 0.1e1 / pkin(1);
t893 = 0.1e1 / pkin(2);
t856 = t882 * t889 + t885 * t892;
t855 = t882 * t892 - t885 * t889;
t854 = t881 * t888 + t884 * t891;
t853 = t881 * t891 - t884 * t888;
t852 = t880 * t887 + t883 * t890;
t851 = t880 * t890 - t883 * t887;
t850 = g(1) * t873 + g(2) * t872;
t849 = g(1) * t872 - g(2) * t873;
t821 = t845 * t882 + t848 * t885;
t820 = -t845 * t885 + t848 * t882;
t819 = t844 * t881 + t847 * t884;
t818 = -t844 * t884 + t847 * t881;
t817 = t843 * t880 + t846 * t883;
t816 = -t843 * t883 + t846 * t880;
t803 = pkin(1) * (-t868 * t882 + t871 * t885) - t827 * pkin(2);
t802 = pkin(1) * (-t867 * t881 + t870 * t884) - t825 * pkin(2);
t801 = pkin(1) * (-t866 * t880 + t869 * t883) - t823 * pkin(2);
t800 = pkin(1) * (t868 * t885 + t871 * t882) + t826 * pkin(2);
t799 = pkin(1) * (t867 * t884 + t870 * t881) + t824 * pkin(2);
t798 = pkin(1) * (t866 * t883 + t869 * t880) + t822 * pkin(2);
t794 = pkin(1) * ((t855 * t873 - t856 * t872) * t871 + t868 * (t855 * t872 + t856 * t873)) - t797 * pkin(2);
t793 = pkin(1) * ((t853 * t873 - t854 * t872) * t870 + t867 * (t853 * t872 + t854 * t873)) - t796 * pkin(2);
t792 = pkin(1) * ((t851 * t873 - t852 * t872) * t869 + t866 * (t851 * t872 + t852 * t873)) - t795 * pkin(2);
t1 = [0, (-t816 * t899 - t818 * t897 - t820 * t895) * t894, (-t817 * t899 - t819 * t897 - t821 * t895) * t894, 0, (-t810 * t899 - t812 * t897 - t814 * t895 + (-t801 * t906 - t802 * t904 - t803 * t902) * t893) * t894, (-t811 * t899 - t813 * t897 - t815 * t895 + (-t801 * t905 - t802 * t903 - t803 * t901) * t893) * t894, 0, 0, 0, -t849 * t872 - t850 * t873; 0, (t816 * t900 + t818 * t898 + t820 * t896) * t894, (t817 * t900 + t819 * t898 + t821 * t896) * t894, 0, (t810 * t900 + t812 * t898 + t814 * t896 + (-t798 * t906 - t799 * t904 - t800 * t902) * t893) * t894, (t811 * t900 + t813 * t898 + t815 * t896 + (-t798 * t905 - t799 * t903 - t800 * t901) * t893) * t894, 0, 0, 0, t849 * t873 - t850 * t872; 0, (-t816 * t909 - t818 * t908 - t820 * t907) * t894, (-t817 * t909 - t819 * t908 - t821 * t907) * t894, 0, (-t795 * t906 - t796 * t904 - t797 * t902 + (-t792 * t906 - t793 * t904 - t794 * t902) * t893) * t894, (-t795 * t905 - t796 * t903 - t797 * t901 + (-t792 * t905 - t793 * t903 - t794 * t901) * t893) * t894, 0, t849, t850, 0;];
tau_reg  = t1;
