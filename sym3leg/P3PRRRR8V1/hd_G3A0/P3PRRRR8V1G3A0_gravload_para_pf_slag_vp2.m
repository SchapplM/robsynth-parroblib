% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G3A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:14
% EndTime: 2020-08-06 17:16:15
% DurationCPUTime: 0.86s
% Computational Cost: add. (675->128), mult. (1503->247), div. (54->7), fcn. (1311->22), ass. (0->108)
t876 = sin(qJ(2,1));
t882 = cos(qJ(2,1));
t870 = legFrame(1,2);
t855 = sin(t870);
t858 = cos(t870);
t839 = g(1) * t858 - g(2) * t855;
t863 = sin(pkin(6));
t865 = cos(pkin(6));
t940 = t865 * g(3);
t833 = -t839 * t863 - t940;
t836 = g(1) * t855 + g(2) * t858;
t864 = sin(pkin(3));
t866 = cos(pkin(3));
t896 = t833 * t866 + t836 * t864;
t941 = g(3) * t863;
t899 = t839 * t865 - t941;
t884 = t896 * t876 + t899 * t882;
t874 = sin(qJ(2,2));
t880 = cos(qJ(2,2));
t869 = legFrame(2,2);
t854 = sin(t869);
t857 = cos(t869);
t838 = g(1) * t857 - g(2) * t854;
t832 = -t838 * t863 - t940;
t835 = g(1) * t854 + g(2) * t857;
t897 = t832 * t866 + t835 * t864;
t900 = t838 * t865 - t941;
t885 = t897 * t874 + t900 * t880;
t872 = sin(qJ(2,3));
t878 = cos(qJ(2,3));
t868 = legFrame(3,2);
t853 = sin(t868);
t856 = cos(t868);
t837 = g(1) * t856 - g(2) * t853;
t831 = -t837 * t863 - t940;
t834 = g(1) * t853 + g(2) * t856;
t898 = t831 * t866 + t834 * t864;
t901 = t837 * t865 - t941;
t886 = t898 * t872 + t901 * t878;
t877 = cos(qJ(3,3));
t944 = pkin(2) * t877;
t840 = -pkin(5) * t878 + t872 * t944;
t871 = sin(qJ(3,3));
t947 = pkin(2) * t871;
t828 = t840 * t864 + t866 * t947;
t950 = 0.1e1 / t828;
t879 = cos(qJ(3,2));
t943 = pkin(2) * t879;
t841 = -pkin(5) * t880 + t874 * t943;
t873 = sin(qJ(3,2));
t946 = pkin(2) * t873;
t829 = t841 * t864 + t866 * t946;
t949 = 0.1e1 / t829;
t881 = cos(qJ(3,1));
t942 = pkin(2) * t881;
t842 = -pkin(5) * t882 + t876 * t942;
t875 = sin(qJ(3,1));
t945 = pkin(2) * t875;
t830 = t842 * t864 + t866 * t945;
t948 = 0.1e1 / t830;
t939 = t950 / t877;
t938 = t949 / t879;
t937 = t948 / t881;
t936 = t950 * t834;
t935 = t949 * t835;
t934 = t948 * t836;
t929 = t834 * t866;
t927 = t835 * t866;
t925 = t836 * t866;
t924 = mrSges(3,2) * t864 * t940;
t923 = t863 * t864;
t922 = t866 * t872;
t921 = t866 * t874;
t920 = t866 * t876;
t919 = t866 * t878;
t918 = t866 * t880;
t917 = t866 * t882;
t916 = t871 * t878;
t915 = t873 * t880;
t914 = t875 * t882;
t913 = (((t831 * t864 - t929) * mrSges(3,1) + t886 * mrSges(3,2)) * t877 + t871 * (t924 + (t837 * t923 + t929) * mrSges(3,2) + t886 * mrSges(3,1))) * t939;
t912 = (((t832 * t864 - t927) * mrSges(3,1) + t885 * mrSges(3,2)) * t879 + t873 * (t924 + (t838 * t923 + t927) * mrSges(3,2) + t885 * mrSges(3,1))) * t938;
t911 = (((t833 * t864 - t925) * mrSges(3,1) + t884 * mrSges(3,2)) * t881 + t875 * (t924 + (t839 * t923 + t925) * mrSges(3,2) + t884 * mrSges(3,1))) * t937;
t867 = mrSges(2,2) - mrSges(3,3);
t910 = (t886 * t867 + (t872 * t901 - t878 * t898) * (mrSges(3,1) * t877 - mrSges(3,2) * t871 + mrSges(2,1))) * t939;
t909 = (t885 * t867 + (t874 * t900 - t880 * t897) * (mrSges(3,1) * t879 - mrSges(3,2) * t873 + mrSges(2,1))) * t938;
t908 = (t884 * t867 + (t876 * t899 - t882 * t896) * (mrSges(3,1) * t881 - mrSges(3,2) * t875 + mrSges(2,1))) * t937;
t907 = ((-t863 * t872 + t865 * t919) * t944 + pkin(5) * (t863 * t878 + t865 * t922)) * t913;
t906 = ((-t863 * t874 + t865 * t918) * t943 + pkin(5) * (t863 * t880 + t865 * t921)) * t912;
t905 = ((-t863 * t876 + t865 * t917) * t942 + pkin(5) * (t863 * t882 + t865 * t920)) * t911;
t892 = t864 * t877 + t871 * t922;
t904 = (t863 * t916 + t865 * t892) * t910;
t891 = t864 * t879 + t873 * t921;
t903 = (t863 * t915 + t865 * t891) * t909;
t890 = t864 * t881 + t875 * t920;
t902 = (t863 * t914 + t865 * t890) * t908;
t895 = -t840 * t866 + t864 * t947;
t894 = -t841 * t866 + t864 * t946;
t893 = -t842 * t866 + t864 * t945;
t883 = 0.1e1 / pkin(2);
t859 = m(1) + m(2) + m(3);
t845 = pkin(5) * t876 + t882 * t942;
t844 = pkin(5) * t874 + t880 * t943;
t843 = pkin(5) * t872 + t878 * t944;
t815 = t865 * t845 + t863 * t893;
t814 = t865 * t844 + t863 * t894;
t813 = t865 * t843 + t863 * t895;
t1 = [-t856 * t904 - t857 * t903 - t858 * t902 - g(1) * m(4) + (-t856 * t907 - t857 * t906 - t858 * t905) * t883 + (-(t815 * t858 + t830 * t855) * t934 - (t814 * t857 + t829 * t854) * t935 - (t813 * t856 + t828 * t853) * t936) * t859; t853 * t904 + t854 * t903 + t855 * t902 - g(2) * m(4) + (t853 * t907 + t854 * t906 + t855 * t905) * t883 + (-(-t815 * t855 + t830 * t858) * t934 - (-t814 * t854 + t829 * t857) * t935 - (-t813 * t853 + t828 * t856) * t936) * t859; (t863 * t890 - t865 * t914) * t908 + (t863 * t891 - t865 * t915) * t909 + (t863 * t892 - t865 * t916) * t910 - g(3) * m(4) + (((t863 * t917 + t865 * t876) * t942 + (t863 * t920 - t865 * t882) * pkin(5)) * t911 + ((t863 * t918 + t865 * t874) * t943 + (t863 * t921 - t865 * t880) * pkin(5)) * t912 + ((t863 * t919 + t865 * t872) * t944 + (t863 * t922 - t865 * t878) * pkin(5)) * t913) * t883 + (-(-t863 * t845 + t865 * t893) * t934 - (-t863 * t844 + t865 * t894) * t935 - (-t863 * t843 + t865 * t895) * t936) * t859;];
taugX  = t1;
