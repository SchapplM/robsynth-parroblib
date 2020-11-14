% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:38
% EndTime: 2020-08-06 18:21:38
% DurationCPUTime: 0.54s
% Computational Cost: add. (376->100), mult. (621->188), div. (63->7), fcn. (654->18), ass. (0->88)
t879 = legFrame(3,3);
t866 = sin(t879);
t869 = cos(t879);
t848 = t866 * g(1) - t869 * g(2);
t851 = t869 * g(1) + t866 * g(2);
t883 = sin(qJ(1,3));
t889 = cos(qJ(1,3));
t821 = t848 * t883 - t851 * t889;
t880 = legFrame(2,3);
t867 = sin(t880);
t870 = cos(t880);
t849 = t867 * g(1) - t870 * g(2);
t852 = t870 * g(1) + t867 * g(2);
t885 = sin(qJ(1,2));
t891 = cos(qJ(1,2));
t823 = t849 * t885 - t852 * t891;
t881 = legFrame(1,3);
t868 = sin(t881);
t871 = cos(t881);
t850 = t868 * g(1) - t871 * g(2);
t853 = t871 * g(1) + t868 * g(2);
t887 = sin(qJ(1,1));
t893 = cos(qJ(1,1));
t825 = t850 * t887 - t853 * t893;
t847 = t868 * t893 + t871 * t887;
t886 = sin(qJ(3,1));
t865 = -t886 * pkin(3) - qJ(2,1);
t862 = 0.1e1 / t865;
t904 = t847 * t862;
t846 = t867 * t891 + t870 * t885;
t884 = sin(qJ(3,2));
t864 = t884 * pkin(3) + qJ(2,2);
t861 = 0.1e1 / t864;
t905 = t846 * t861;
t845 = t866 * t889 + t869 * t883;
t882 = sin(qJ(3,3));
t863 = t882 * pkin(3) + qJ(2,3);
t860 = 0.1e1 / t863;
t906 = t845 * t860;
t916 = t821 * t906 + t823 * t905 - t825 * t904;
t844 = -t868 * t887 + t871 * t893;
t907 = t844 * t862;
t843 = -t867 * t885 + t870 * t891;
t908 = t843 * t861;
t842 = -t866 * t883 + t869 * t889;
t909 = t842 * t860;
t915 = t821 * t909 + t823 * t908 - t825 * t907;
t822 = t848 * t889 + t851 * t883;
t824 = t849 * t891 + t852 * t885;
t826 = t850 * t893 + t853 * t887;
t914 = t822 * t906 + t824 * t905 - t826 * t904;
t913 = t822 * t909 + t824 * t908 - t826 * t907;
t854 = -t883 * g(1) + t889 * g(2);
t855 = t889 * g(1) + t883 * g(2);
t912 = (-t854 * t866 - t855 * t869) * t860;
t856 = -t885 * g(1) + t891 * g(2);
t857 = t891 * g(1) + t885 * g(2);
t911 = (-t856 * t867 - t857 * t870) * t861;
t858 = g(1) * t887 - g(2) * t893;
t859 = g(1) * t893 + g(2) * t887;
t910 = (t858 * t868 - t859 * t871) * t862;
t903 = t882 * t912;
t888 = cos(qJ(3,3));
t902 = t888 * t912;
t901 = t884 * t911;
t890 = cos(qJ(3,2));
t900 = t890 * t911;
t899 = t886 * t910;
t892 = cos(qJ(3,1));
t898 = t892 * t910;
t897 = t854 * t869 - t855 * t866;
t896 = t856 * t870 - t857 * t867;
t895 = t858 * t871 + t859 * t868;
t894 = 0.1e1 / pkin(3);
t878 = 0.1e1 / t886;
t877 = 0.1e1 / t884;
t876 = 0.1e1 / t882;
t875 = pkin(1) + pkin(5) + pkin(6);
t841 = t865 * t893 + t875 * t887;
t840 = -t864 * t891 + t875 * t885;
t839 = -t863 * t889 + t875 * t883;
t838 = -t887 * t865 + t875 * t893;
t837 = t885 * t864 + t875 * t891;
t836 = t883 * t863 + t875 * t889;
t820 = -t853 * (-t887 * pkin(1) + t893 * qJ(2,1)) + t850 * (t893 * pkin(1) + t887 * qJ(2,1));
t819 = -t852 * (-t885 * pkin(1) + t891 * qJ(2,2)) + t849 * (t891 * pkin(1) + t885 * qJ(2,2));
t818 = -t851 * (-t883 * pkin(1) + t889 * qJ(2,3)) + t848 * (t889 * pkin(1) + t883 * qJ(2,3));
t1 = [0, t913, -t915, -t913, t915, -(t844 * t820 - (t838 * t871 - t868 * t841) * t826) * t862 + (t843 * t819 - (t837 * t870 - t867 * t840) * t824) * t861 + (t842 * t818 - (t836 * t869 - t866 * t839) * t822) * t860, 0, 0, 0, 0, 0, t842 * t903 + t843 * t901 - t844 * t899, t842 * t902 + t843 * t900 - t844 * t898, -g(1); 0, t914, -t916, -t914, t916, -(t847 * t820 - (t868 * t838 + t841 * t871) * t826) * t862 + (t846 * t819 - (t867 * t837 + t840 * t870) * t824) * t861 + (t845 * t818 - (t866 * t836 + t839 * t869) * t822) * t860, 0, 0, 0, 0, 0, t845 * t903 + t846 * t901 - t847 * t899, t845 * t902 + t846 * t900 - t847 * t898, -g(2); 0, 0, 0, 0, 0, -t822 * t876 * t888 - t824 * t877 * t890 - t826 * t878 * t892, 0, 0, 0, 0, 0, (-t878 * (g(3) * t886 - t892 * t895) - t877 * (g(3) * t884 + t890 * t896) - t876 * (g(3) * t882 + t888 * t897)) * t894, (-t878 * (g(3) * t892 + t886 * t895) - t877 * (g(3) * t890 - t884 * t896) - t876 * (g(3) * t888 - t882 * t897)) * t894, -g(3);];
tau_reg  = t1;
