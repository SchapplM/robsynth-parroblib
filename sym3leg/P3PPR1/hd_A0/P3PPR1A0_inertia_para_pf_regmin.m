% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PPR1A0
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
% tau_reg [3*(3+1)/2x6]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PPR1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_inertia_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_inertia_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (150->50), mult. (268->77), div. (0->0), fcn. (266->8), ass. (0->44)
t108 = legFrame(3,3);
t100 = sin(t108);
t109 = legFrame(2,3);
t101 = sin(t109);
t110 = legFrame(1,3);
t102 = sin(t110);
t103 = cos(t108);
t104 = cos(t109);
t105 = cos(t110);
t111 = xP(3);
t106 = sin(t111);
t107 = cos(t111);
t112 = koppelP(3,2);
t115 = koppelP(3,1);
t81 = t100 * t115 - t103 * t112;
t82 = t100 * t112 + t103 * t115;
t74 = t106 * t81 + t82 * t107;
t113 = koppelP(2,2);
t116 = koppelP(2,1);
t83 = t101 * t116 - t104 * t113;
t84 = t101 * t113 + t104 * t116;
t76 = t106 * t83 + t84 * t107;
t114 = koppelP(1,2);
t117 = koppelP(1,1);
t85 = t102 * t117 - t105 * t114;
t86 = t102 * t114 + t105 * t117;
t78 = t106 * t85 + t86 * t107;
t87 = -t106 * t115 - t107 * t112;
t88 = -t106 * t116 - t107 * t113;
t89 = -t106 * t117 - t107 * t114;
t90 = -t106 * t112 + t107 * t115;
t91 = -t106 * t113 + t107 * t116;
t92 = -t106 * t114 + t107 * t117;
t122 = t78 * (-t89 * t102 + t92 * t105) + t76 * (-t88 * t101 + t91 * t104) + t74 * (-t87 * t100 + t90 * t103);
t121 = t74 * t103 + t76 * t104 + t78 * t105;
t120 = t100 ^ 2 + t101 ^ 2 + t102 ^ 2;
t119 = t103 ^ 2 + t104 ^ 2 + t105 ^ 2;
t118 = -t74 * t100 - t76 * t101 - t78 * t102;
t93 = t106 ^ 2 + t107 ^ 2;
t80 = t119 + t120;
t79 = -t86 * t106 + t85 * t107;
t77 = -t84 * t106 + t83 * t107;
t75 = -t82 * t106 + t81 * t107;
t1 = [t120, t80, 0, 0, 0, t93; -t103 * t100 - t104 * t101 - t105 * t102, 0, 0, 0, 0, 0; t119, t80, 0, 0, 0, t93; t118, t75 * t103 + t77 * t104 + t79 * t105 + t118, 0, -t106, -t107, 0; t121, t75 * t100 + t77 * t101 + t79 * t102 + t121, 0, t107, -t106, 0; t122, t79 * (t92 * t102 + t89 * t105) + t77 * (t91 * t101 + t88 * t104) + t75 * (t90 * t100 + t87 * t103) + t122, 1, 0, 0, 0;];
tau_reg  = t1;
