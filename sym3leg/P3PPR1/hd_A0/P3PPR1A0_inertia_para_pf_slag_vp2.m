% Calculate inertia matrix for parallel robot
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PPR1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:30
% EndTime: 2019-05-03 14:37:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (173->44), mult. (296->66), div. (0->0), fcn. (248->8), ass. (0->39)
t91 = legFrame(3,3);
t83 = sin(t91);
t92 = legFrame(2,3);
t84 = sin(t92);
t93 = legFrame(1,3);
t85 = sin(t93);
t106 = t83 ^ 2 + t84 ^ 2 + t85 ^ 2;
t86 = cos(t91);
t87 = cos(t92);
t88 = cos(t93);
t105 = t86 ^ 2 + t87 ^ 2 + t88 ^ 2;
t103 = koppelP(1,1);
t102 = koppelP(2,1);
t101 = koppelP(3,1);
t100 = koppelP(1,2);
t99 = koppelP(2,2);
t98 = koppelP(3,2);
t97 = mrSges(3,1);
t96 = mrSges(3,2);
t95 = xP(3);
t94 = m(1) + m(2);
t90 = cos(t95);
t89 = sin(t95);
t76 = -t100 * t89 + t103 * t90;
t75 = t102 * t90 - t89 * t99;
t74 = t101 * t90 - t89 * t98;
t73 = -t100 * t90 - t103 * t89;
t72 = -t102 * t89 - t90 * t99;
t71 = -t101 * t89 - t90 * t98;
t70 = t73 * t88 + t76 * t85;
t69 = -t73 * t85 + t76 * t88;
t68 = t72 * t87 + t75 * t84;
t67 = -t72 * t84 + t75 * t87;
t66 = t71 * t86 + t74 * t83;
t65 = -t71 * t83 + t74 * t86;
t64 = (m(2) - t94) * (t83 * t86 + t84 * t87 + t85 * t88);
t63 = -t89 * t96 + t90 * t97 + (t65 * t86 + t67 * t87 + t69 * t88) * t94 + (t66 * t83 + t68 * t84 + t70 * t85) * m(2);
t62 = -t89 * t97 - t90 * t96 + (-t65 * t83 - t67 * t84 - t69 * t85) * t94 + (t66 * t86 + t68 * t87 + t70 * t88) * m(2);
t1 = [t105 * m(2) + t106 * t94 + m(3), t64, t62; t64, t106 * m(2) + t105 * t94 + m(3), t63; t62, t63, Ifges(3,3) + (t65 ^ 2 + t67 ^ 2 + t69 ^ 2) * t94 + (t66 ^ 2 + t68 ^ 2 + t70 ^ 2) * m(2);];
MX  = t1;
