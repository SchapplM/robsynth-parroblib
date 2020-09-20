% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PPR1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x6]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PPR1G1P1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1P1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (132->42), mult. (359->75), div. (0->0), fcn. (238->8), ass. (0->43)
t159 = xP(3);
t153 = sin(t159);
t155 = xDP(3) ^ 2;
t170 = t153 * t155;
t154 = cos(t159);
t169 = t154 * t155;
t160 = koppelP(3,2);
t163 = koppelP(3,1);
t141 = t153 * t163 + t154 * t160;
t144 = -t153 * t160 + t154 * t163;
t156 = legFrame(3,3);
t147 = sin(t156);
t150 = cos(t156);
t129 = (-t141 * t150 + t144 * t147) * t155;
t161 = koppelP(2,2);
t164 = koppelP(2,1);
t142 = t153 * t164 + t154 * t161;
t145 = -t153 * t161 + t154 * t164;
t157 = legFrame(2,3);
t148 = sin(t157);
t151 = cos(t157);
t130 = (-t142 * t151 + t145 * t148) * t155;
t162 = koppelP(1,2);
t165 = koppelP(1,1);
t143 = t153 * t165 + t154 * t162;
t146 = -t153 * t162 + t154 * t165;
t158 = legFrame(1,3);
t149 = sin(t158);
t152 = cos(t158);
t131 = (-t143 * t152 + t146 * t149) * t155;
t135 = t147 * t163 - t150 * t160;
t136 = t147 * t160 + t150 * t163;
t137 = t148 * t164 - t151 * t161;
t138 = t148 * t161 + t151 * t164;
t139 = t149 * t165 - t152 * t162;
t140 = t149 * t162 + t152 * t165;
t168 = (t153 * t139 + t140 * t154) * t131 + (t153 * t137 + t138 * t154) * t130 + (t153 * t135 + t136 * t154) * t129;
t167 = t150 * t129 + t151 * t130 + t152 * t131;
t166 = -t147 * t129 - t148 * t130 - t149 * t131;
t134 = (-t143 * t149 - t146 * t152) * t155;
t133 = (-t142 * t148 - t145 * t151) * t155;
t132 = (-t141 * t147 - t144 * t150) * t155;
t1 = [t166, t150 * t132 + t151 * t133 + t152 * t134 + t166, 0, -t169, t170, 0; t167, t147 * t132 + t148 * t133 + t149 * t134 + t167, 0, -t170, -t169, 0; t168, (t139 * t154 - t140 * t153) * t134 + (t137 * t154 - t138 * t153) * t133 + (t135 * t154 - t136 * t153) * t132 + t168, 0, 0, 0, 0;];
tau_reg  = t1;
