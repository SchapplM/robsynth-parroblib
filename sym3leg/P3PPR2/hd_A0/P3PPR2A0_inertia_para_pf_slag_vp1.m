% Calculate inertia matrix for parallel robot
% P3PPR2A0
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
%   mass of all robot links (including platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:31
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MX = P3PPR2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:31:32
% EndTime: 2018-12-20 17:31:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (175->46), mult. (303->71), div. (0->0), fcn. (248->8), ass. (0->39)
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
t97 = rSges(3,1);
t96 = rSges(3,2);
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
t69 = t73 * t85 - t76 * t88;
t68 = t72 * t87 + t75 * t84;
t67 = t72 * t84 - t75 * t87;
t66 = t71 * t86 + t74 * t83;
t65 = t71 * t83 - t74 * t86;
t64 = (-m(2) + t94) * (t83 * t86 + t84 * t87 + t85 * t88);
t63 = m(3) * (-t89 * t97 - t90 * t96) + (t86 * t66 + t87 * t68 + t88 * t70) * t94 + (t83 * t65 + t84 * t67 + t85 * t69) * m(2);
t62 = m(3) * (-t89 * t96 + t90 * t97) + (t83 * t66 + t84 * t68 + t85 * t70) * t94 + (-t86 * t65 - t87 * t67 - t88 * t69) * m(2);
t1 = [t106 * m(2) + t105 * t94 + m(3), t64, t63; t64, t105 * m(2) + t106 * t94 + m(3), t62; t63, t62, Icges(3,3) + m(3) * (t96 ^ 2 + t97 ^ 2) + (t66 ^ 2 + t68 ^ 2 + t70 ^ 2) * t94 + (t65 ^ 2 + t67 ^ 2 + t69 ^ 2) * m(2);];
MX  = t1;
