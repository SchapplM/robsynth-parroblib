% Calculate vector of centrifugal and coriolis load on the joints for
% P3PPR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PPR1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:30
% EndTime: 2019-05-03 14:37:31
% DurationCPUTime: 0.17s
% Computational Cost: add. (109->45), mult. (233->89), div. (0->0), fcn. (160->8), ass. (0->46)
t140 = xDP(3);
t136 = t140 ^ 2;
t142 = xP(3);
t134 = sin(t142);
t135 = cos(t142);
t145 = koppelP(3,2);
t148 = koppelP(3,1);
t116 = t134 * t148 + t135 * t145;
t119 = -t134 * t145 + t135 * t148;
t137 = legFrame(3,3);
t128 = sin(t137);
t131 = cos(t137);
t157 = m(2) * (-t116 * t131 + t119 * t128);
t146 = koppelP(2,2);
t149 = koppelP(2,1);
t117 = t134 * t149 + t135 * t146;
t120 = -t134 * t146 + t135 * t149;
t138 = legFrame(2,3);
t129 = sin(t138);
t132 = cos(t138);
t156 = m(2) * (-t117 * t132 + t120 * t129);
t147 = koppelP(1,2);
t150 = koppelP(1,1);
t118 = t134 * t150 + t135 * t147;
t121 = -t134 * t147 + t135 * t150;
t139 = legFrame(1,3);
t130 = sin(t139);
t133 = cos(t139);
t155 = m(2) * (-t118 * t133 + t121 * t130);
t141 = m(1) + m(2);
t154 = m(2) - t141;
t153 = (t116 * t128 + t119 * t131) * t141;
t152 = (t117 * t129 + t120 * t132) * t141;
t151 = (t118 * t130 + t121 * t133) * t141;
t144 = mrSges(3,1);
t143 = mrSges(3,2);
t127 = t133 ^ 2;
t126 = t132 ^ 2;
t125 = t131 ^ 2;
t124 = t130 ^ 2;
t123 = t129 ^ 2;
t122 = t128 ^ 2;
t115 = t154 * t133 * t130;
t114 = t154 * t132 * t129;
t113 = t154 * t131 * t128;
t1 = [(-(t127 * m(2) + t124 * t141) * t121 - t115 * t118 - (t126 * m(2) + t123 * t141) * t120 - t114 * t117 - (t125 * m(2) + t122 * t141) * t119 - t113 * t116 + t134 * t143 - t135 * t144) * t136; (-t115 * t121 - (t124 * m(2) + t127 * t141) * t118 - t114 * t120 - (t123 * m(2) + t126 * t141) * t117 - t113 * t119 - (t122 * m(2) + t125 * t141) * t116 - t134 * t144 - t135 * t143) * t136; (-(-t130 * t151 + t133 * t155) * t121 - (t130 * t155 + t133 * t151) * t118 - (-t129 * t152 + t132 * t156) * t120 - (t129 * t156 + t132 * t152) * t117 - (-t128 * t153 + t131 * t157) * t119 - (t128 * t157 + t131 * t153) * t116) * t136;];
taucX  = t1;
