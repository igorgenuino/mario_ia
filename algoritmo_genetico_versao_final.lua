--- Projeto Inteligencia Artificial na Robotica
--- 
---Igor Genuino
---
--- 
--- versão do codigo: 16

Pasta = "fase01" ----- Folder caso desejamos aplicar esse algoritmo em outro jogos

if gameinfo.getromname() == "Super Mario World (USA)" then
	Diretorio = Pasta
	Estado = "estado_inicial" ---- estado inicial / pode ser mudado por fase se desejar
	ButtonNames = { --nome dos botões conforme controle do nintendo
		"A",
		"B",
		"X",
		"Y",
		"Up",
		"Down",
		"Left",
		"Right",
	}
elseif gameinfo.getromname() == "Super Mario Bros." then --- caso sua versão da row não seja americana
	Diretorio = "SMB1"
	Estado = "SMB1-1" --estado do jogo salvo
	ButtonNames = {
		"A",
		"B",
		"Up",
		"Down",
		"Left",
		"Right",
	}
end

Diretorio_Armazenamento = Pasta ---- define como diretorio para arquivar os dados 

BoxRadius = 6 ------ funcao que mapea os objetos na fase
InputSize = (BoxRadius*2+1)*(BoxRadius*2+1)  ----- gera o comprimento do BoxRadius "nao entendi ainda rs "

Inputs = InputSize+1 ------- Leva esse comprimento para exibir no quadro
Outputs = #ButtonNames ----- representa a entrada do botões controlado pelo algoritmo

Populacao = 300  ----- população usada a cada geração
D_similaridade = 2.0 ----- variacao da Similaridade entras as especies , elemento disjuntos
D_pesos = 0.4 ------ variacao do peso da rede neural
D_limiar = 1.0 ------ variacao de do limiar 

Especies_Salvas = 15 ----- quantidade de especies salvas após fim de uma geracao

Mutacoes_conectadas = 0.25 -----  taxa de mutações relacionadas aos indiviuos 
Pertubacao = 0.90 ----- 
Especies_mortas = 0.75 ----- taxa de quantas especies são eliminadas
Combinacao_mutacoes = 2.0 ----- taxa com que combinamos as especies
variacao_camadas = 0.50 ----- taxa de vari
variacao_camada_escondida = 0.40 -----
taxa_apreendizagem = 0.1 ------ defini velocidade que a rede aprende
Desabilita_mutacoes = 0.4 ------ taxa de chance de se ter mutações desabilitadas, importante por que nosso cenário muda
Habilita_mutacoes = 0.2 -------

Tempo = 20 -------- defini um tempo para que se nada ocorra, recomeça o codigo

MaxNodes = 1000000 -------- numero de pontos maximo para o reconhecimento de padrão da rede

estado_pool = {}
clock = os.clock()
frameCount = emu.framecount()

function getPositions() ------ configuraçao da posição do mario durante a jogo, funçao do emulador
	if gameinfo.getromname() == "Super Mario World (USA)" then
		marioX = memory.read_s16_le(0x94) 
		marioY = memory.read_s16_le(0x96)

		local layer1x = memory.read_s16_le(0x1A)
		local layer1y = memory.read_s16_le(0x1C)

		screenX = marioX-layer1x
		screenY = marioY-layer1y
	elseif gameinfo.getromname() == "Super Mario Bros." then
		marioX = memory.readbyte(0x6D) * 0x100 + memory.readbyte(0x86)
		marioY = memory.readbyte(0x03B8)+16

		screenX = memory.readbyte(0x03AD)
		screenY = memory.readbyte(0x03B8)
	end
end

function getTile(dx, dy) ------ função para mapear o cenario, função emulador
	if gameinfo.getromname() == "Super Mario World (USA)" then
		x = math.floor((marioX+dx+8)/16)
		y = math.floor((marioY+dy)/16)

		return memory.readbyte(0x1C800 + math.floor(x/0x10)*0x1B0 + y*0x10 + x%0x10)
	
	
	
	
	elseif gameinfo.getromname() == "Super Mario Bros." then
		local x = marioX + dx + 8
		local y = marioY + dy - 16
		local page = math.floor(x/256)%2

		local subx = math.floor((x%256)/16)
		local suby = math.floor((y - 32)/16)
		local addr = 0x500 + page*13*16+suby*16+subx

		if suby >= 13 or suby < 0 then
			return 0
		end

		if memory.readbyte(addr) ~= 0 then
			return 1
		else
			return 0
		end
	end
end

function getSprites() ---função para mapear os inimigos no jogo e cenario, função do emulador
	if gameinfo.getromname() == "Super Mario World (USA)" then
		local sprites = {}
		for slot=0,11 do
			local status = memory.readbyte(0x14C8+slot)
			if status ~= 0 then
				spritex = memory.readbyte(0xE4+slot) + memory.readbyte(0x14E0+slot)*256
				spritey = memory.readbyte(0xD8+slot) + memory.readbyte(0x14D4+slot)*256
				sprites[#sprites+1] = {["x"]=spritex, ["y"]=spritey}
			end
		end

		return sprites
	elseif gameinfo.getromname() == "Super Mario Bros." then
		local sprites = {}
		for slot=0,4 do
			local enemy = memory.readbyte(0xF+slot)
			if enemy ~= 0 then
				local ex = memory.readbyte(0x6E + slot)*0x100 + memory.readbyte(0x87+slot)
				local ey = memory.readbyte(0xCF + slot)+24
				sprites[#sprites+1] = {["x"]=ex,["y"]=ey}
			end
		end

		return sprites
	end
end


function getExtendedSprites() ---função para mapear a posição extendida, conforme se movimentão, dos inimigos, função do emulador
	if gameinfo.getromname() == "Super Mario World (USA)" then
		local extended = {}
		for slot=0,11 do
			local number = memory.readbyte(0x170B+slot)
			if number ~= 0 then
				spritex = memory.readbyte(0x171F+slot) + memory.readbyte(0x1733+slot)*256
				spritey = memory.readbyte(0x1715+slot) + memory.readbyte(0x1729+slot)*256
				extended[#extended+1] = {["x"]=spritex, ["y"]=spritey}
			end
		end

		return extended
	elseif gameinfo.getromname() == "Super Mario Bros." then
		return {}
	end
end

function getInputs() --função do emulador que recebe as entradas do jogo, recebe o mapeamento do cenario, inimigos e mario, 
	getPositions()

	sprites = getSprites() 
	extended = getExtendedSprites()

	local inputs = {}

	for dy=-BoxRadius*16,BoxRadius*16,16 do
		for dx=-BoxRadius*16,BoxRadius*16,16 do
			inputs[#inputs+1] = 0

			tile = getTile(dx, dy)
			if tile == 1 and marioY+dy < 0x1B0 then
				inputs[#inputs] = 1
			end

			for i = 1,#sprites do
				distx = math.abs(sprites[i]["x"] - (marioX+dx))
				disty = math.abs(sprites[i]["y"] - (marioY+dy))
				if distx <= 8 and disty <= 8 then
					inputs[#inputs] = -1
				end
			end

			for i = 1,#extended do
				distx = math.abs(extended[i]["x"] - (marioX+dx))
				disty = math.abs(extended[i]["y"] - (marioY+dy))
				if distx < 8 and disty < 8 then
					inputs[#inputs] = -1
				end
			end
		end
	end

	--mariovx = memory.read_s8(0x7B)
	--mariovy = memory.read_s8(0x7D)

	return inputs
end

function isDead() -- função do emulador que determina "fim da vida do mario"
	if gameinfo.getromname() == "Super Mario World (USA)" then
		return memory.readbyte(0x13e0) == 0x3e
	end
end

------------------------------------------------ Algoritmo Genético Hibrido -------------------------------------------------------------



function Sigmoidal(x) ---- função Logistica Sigmoidal
	return 2/(1+math.exp(-4.9*x))-1
end

function Funcao_distancia() -----  função que definir a distancia percorrida pelo mario
	return rightmost - initialRightMost
end

function Novo() ------ faz um incremento da variavel pool.novo afim gera um nova geracao
	pool.novo = pool.novo + 1
	return pool.novo
end

function novaPool()  ---função que salva as informações e individuos do algoritmo genetico
	local pool = {}
	pool.especies = {}
	pool.geracao = 0
	pool.novo = Outputs
	pool.especie_atual = 1
	pool.genoma_atual = 1
	pool.quadro_atual = 0
	pool.maxima_distancia = 0

	return pool
end

function novas_especies() ----função que avalia a cada iteração todos individuos e registra os indivuos conforme abaixo
	local especies = {}
	especies.melhor_distancia = 0
	especies.velha = 0
	especies.genomes = {}
	especies.Media_distancia = 0

	return especies
end

function newGenome() ---função que atribui os novos valores ao genoma usando as taxas definidas 
	local genome = {}
	genome.genes = {}
	genome.distancia = 0
	genome.adjusteddistancia = 0
	genome.network = {}
	genome.maxneuron = 0
	genome.globalRank = 0
	genome.mutacoes_valores = {}
	genome.mutacoes_valores["conexoes"] = Mutacoes_conectadas
	genome.mutacoes_valores["ligacao"] = Combinacao_mutacoes
	genome.mutacoes_valores["tendencia"] = variacao_camada_escondida
	genome.mutacoes_valores["no"] = variacao_camadas
	genome.mutacoes_valores["habilita"] = Habilita_mutacoes
	genome.mutacoes_valores["incapacitar"] = Desabilita_mutacoes
	genome.mutacoes_valores["passo"] = taxa_apreendizagem

	return genome
end

function copiaGenome(genome) ---essa função faz uma copia do genoma com as taxas para comparação 
	local genome2 = newGenome()
	for g=1,#genome.genes do
		table.insert(genome2.genes, copyGene(genome.genes[g]))
	end
	genome2.maxneuron = genome.maxneuron
	genome2.mutacoes_valores["conexoes"] = genome.mutacoes_valores["conexoes"]
	genome2.mutacoes_valores["ligacao"] = genome.mutacoes_valores["ligacao"]
	genome2.mutacoes_valores["tendencia"] = genome.mutacoes_valores["tendencia"]
	genome2.mutacoes_valores["no"] = genome.mutacoes_valores["no"]
	genome2.mutacoes_valores["habilita"] = genome.mutacoes_valores["habilita"]
	genome2.mutacoes_valores["incapacitar"] = genome.mutacoes_valores["incapacitar"]

	return genome2
end

function Genome_basico() -----essa função utiliza o genoma novo com as taxas aplicadas para treinamento da rede
	local genome = newGenome()
	local novo = 1

	genome.maxneuron = Inputs
	mutate(genome)

	return genome
end

function newGene()
	local gene = {}
	gene.into = 0
	gene.out = 0
	gene.weight = 0.0
	gene.enabled = true
	gene.novo = 0

	return gene
end

function copyGene(gene)
	local gene2 = newGene()
	gene2.into = gene.into
	gene2.out = gene.out
	gene2.weight = gene.weight
	gene2.enabled = gene.enabled
	gene2.novo = gene.novo

	return gene2
end

function newNeuron() ---defini um neuronio 
	local neuron = {}
	neuron.incoming = {}
	neuron.value = 0.0

	return neuron
end

function generateNetwork(genome) ---função gera a rede neural, como entrada os individuos genoma do algoritmo genetico
	local network = {}
	network.neurons = {}

	for i=1,Inputs do
		network.neurons[i] = newNeuron()
	end

	for o=1,Outputs do
		network.neurons[MaxNodes+o] = newNeuron()
	end

	table.sort(genome.genes, function (a,b) 
		return (a.out < b.out)
	end)
	for i=1,#genome.genes do
		local gene = genome.genes[i]
		if gene.enabled then
			if network.neurons[gene.out] == nil then
				network.neurons[gene.out] = newNeuron()
			end
			local neuron = network.neurons[gene.out]
			table.insert(neuron.incoming, gene)
			if network.neurons[gene.into] == nil then
				network.neurons[gene.into] = newNeuron()
			end
		end
	end

	genome.network = network
end

function evaluateNetwork(network, inputs) ---função que avalia a rede neural  utiliza a sigmoide para classificar
	table.insert(inputs, 1)
	if #inputs ~= Inputs then
		console.writeline("Incorrect number of neural network inputs.")
		return {}
	end

	for i=1,Inputs do
		network.neurons[i].value = inputs[i]
	end

	for _,neuron in pairs(network.neurons) do
		local sum = 0
		for j = 1,#neuron.incoming do
			local incoming = neuron.incoming[j]
			local other = network.neurons[incoming.into]
			sum = sum + incoming.weight * other.value
		end

		if #neuron.incoming > 0 then
			neuron.value = Sigmoidal(sum)
		end
	end

	local outputs = {}
	for o=1,Outputs do
		local button = "P1 " .. ButtonNames[o]
		if network.neurons[MaxNodes+o].value > 0 then
			outputs[button] = true
		else
			outputs[button] = false
		end
	end

	return outputs
end

function crossover(g1, g2) --- crossover do algoritmo genetico, gera os proximos individuos, filhos
	if g2.distancia > g1.distancia then
		tempg = g1
		g1 = g2
		g2 = tempg
	end

	local child = newGenome()

	local innovations2 = {}
	for i=1,#g2.genes do
		local gene = g2.genes[i]
		innovations2[gene.novo] = gene
	end

	for i=1,#g1.genes do
		local gene1 = g1.genes[i]
		local gene2 = innovations2[gene1.novo]
		if gene2 ~= nil and math.random(2) == 1 and gene2.enabled then
			table.insert(child.genes, copyGene(gene2))
		else
			table.insert(child.genes, copyGene(gene1))
		end
	end

	child.maxneuron = math.max(g1.maxneuron,g2.maxneuron)

	for mutation,rate in pairs(g1.mutacoes_valores) do
		child.mutacoes_valores[mutation] = rate
	end

	return child
end

function randomNeuron(genes, nonInput) -- função que randomiza os neuronios da rede neural para otimizar o algoritmo genetico
	local neurons = {}
	if not nonInput then
		for i=1,Inputs do
			neurons[i] = true
		end
	end
	for o=1,Outputs do
		neurons[MaxNodes+o] = true
	end
	for i=1,#genes do
		if (not nonInput) or genes[i].into > Inputs then
			neurons[genes[i].into] = true
		end
		if (not nonInput) or genes[i].out > Inputs then
			neurons[genes[i].out] = true
		end
	end

	local count = 0
	for _,_ in pairs(neurons) do
		count = count + 1
	end
	local n = math.random(1, count)

	for k,v in pairs(neurons) do
		n = n-1
		if n == 0 then
			return k
		end
	end

	return 0
end

function containsLink(genes, link)
	for i=1,#genes do
		local gene = genes[i]
		if gene.into == link.into and gene.out == link.out then
			return true
		end
	end
end

function pointMutate(genome)
	local step = genome.mutacoes_valores["passo"]

	for i=1,#genome.genes do
		local gene = genome.genes[i]
		if math.random() < Pertubacao then
			gene.weight = gene.weight + math.random() * step*2 - step
		else
			gene.weight = math.random()*4-2
		end
	end
end

function linkMutate(genome, forceBias)
	local neuron1 = randomNeuron(genome.genes, false)
	local neuron2 = randomNeuron(genome.genes, true)

	local newLink = newGene()
	if neuron1 <= Inputs and neuron2 <= Inputs then
		--Both input nodes
		return
	end
	if neuron2 <= Inputs then
		-- Swap output and input
		local temp = neuron1
		neuron1 = neuron2
		neuron2 = temp
	end

	newLink.into = neuron1
	newLink.out = neuron2
	if forceBias then
		newLink.into = Inputs
	end

	if containsLink(genome.genes, newLink) then
		return
	end
	newLink.novo = Novo()
	newLink.weight = math.random()*4-2

	table.insert(genome.genes, newLink)
end

function nodeMutate(genome)  ---faz a mutação pontual 
	if #genome.genes == 0 then
		return
	end

	genome.maxneuron = genome.maxneuron + 1

	local gene = genome.genes[math.random(1,#genome.genes)]
	if not gene.enabled then
		return
	end
	gene.enabled = false

	local gene1 = copyGene(gene)
	gene1.out = genome.maxneuron
	gene1.weight = 1.0
	gene1.novo = Novo()
	gene1.enabled = true
	table.insert(genome.genes, gene1)

	local gene2 = copyGene(gene)
	gene2.into = genome.maxneuron
	gene2.novo = Novo()
	gene2.enabled = true
	table.insert(genome.genes, gene2)
end

function enableDisableMutate(genome, enable) --avalia o potencial de se inserir ou nao mutação no genoma
	local candidates = {}
	for _,gene in pairs(genome.genes) do
		if gene.enabled == not enable then
			table.insert(candidates, gene)
		end
	end

	if #candidates == 0 then
		return
	end

	local gene = candidates[math.random(1,#candidates)]
	gene.enabled = not gene.enabled
end

function mutate(genome) --faz a mutação no genoma
	for mutation,rate in pairs(genome.mutacoes_valores) do
		if math.random(1,2) == 1 then
			genome.mutacoes_valores[mutation] = 0.95*rate
		else
			genome.mutacoes_valores[mutation] = 1.05263*rate
		end
	end

	if math.random() < genome.mutacoes_valores["conexoes"] then
		pointMutate(genome)
	end

	local p = genome.mutacoes_valores["ligacao"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, false)
		end
		p = p - 1
	end

	p = genome.mutacoes_valores["tendencia"]
	while p > 0 do
		if math.random() < p then
			linkMutate(genome, true)
		end
		p = p - 1
	end

	p = genome.mutacoes_valores["no"]
	while p > 0 do
		if math.random() < p then
			nodeMutate(genome)
		end
		p = p - 1
	end

	p = genome.mutacoes_valores["habilita"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, true)
		end
		p = p - 1
	end

	p = genome.mutacoes_valores["incapacitar"]
	while p > 0 do
		if math.random() < p then
			enableDisableMutate(genome, false)
		end
		p = p - 1
	end
end

function disjoint(genes1, genes2)
	local i1 = {}
	for i = 1,#genes1 do
		local gene = genes1[i]
		i1[gene.novo] = true
	end

	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.novo] = true
	end

	local disjointGenes = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if not i2[gene.novo] then
			disjointGenes = disjointGenes+1
		end
	end

	for i = 1,#genes2 do
		local gene = genes2[i]
		if not i1[gene.novo] then
			disjointGenes = disjointGenes+1
		end
	end

	local n = math.max(#genes1, #genes2)

	return disjointGenes / n
end

function weights(genes1, genes2)
	local i2 = {}
	for i = 1,#genes2 do
		local gene = genes2[i]
		i2[gene.novo] = gene
	end

	local sum = 0
	local coincident = 0
	for i = 1,#genes1 do
		local gene = genes1[i]
		if i2[gene.novo] ~= nil then
			local gene2 = i2[gene.novo]
			sum = sum + math.abs(gene.weight - gene2.weight)
			coincident = coincident + 1
		end
	end

	return sum / coincident
end

function sameSpecies(genome1, genome2)
	local dd = D_similaridade*disjoint(genome1.genes, genome2.genes)
	local dw = D_pesos*weights(genome1.genes, genome2.genes)
	return dd + dw < D_limiar
end

function rankGlobally() ---classifica globalmente os individuos
	local global = {}
	for s = 1,#pool.especies do
		local especies = pool.especies[s]
		for g = 1,#especies.genomes do
			table.insert(global, especies.genomes[g])
		end
	end
	table.sort(global, function (a,b)
		return (a.distancia < b.distancia)
	end)

	for g=1,#global do
		global[g].globalRank = g
	end
end

function Calculo_media_distancia(especies) --calcula a media da distancia dos individuos para classificar
	local total = 0

	for g=1,#especies.genomes do
		local genome = especies.genomes[g]
		total = total + genome.globalRank
	end

	especies.Media_distancia = total / #especies.genomes
end

function Distancia_total_media() ---avalia a distancia total percorrida 
	local total = 0
	for s = 1,#pool.especies do
		local especies = pool.especies[s]
		total = total + especies.Media_distancia
	end

	return total
end

function cullSpecies(cutToOne) --corta os individuos muito abaixo do esperado
	for s = 1,#pool.especies do
		local especies = pool.especies[s]

		table.sort(especies.genomes, function (a,b)
			return (a.distancia > b.distancia)
		end)

		local remaining = math.ceil(#especies.genomes/2)
		if cutToOne then
			remaining = 1
		end
		while #especies.genomes > remaining do
			table.remove(especies.genomes)
		end
	end
end

function breedChild(especies)  ----criação de novas especies 
	local child = {}
	if math.random() < Especies_mortas then
		g1 = especies.genomes[math.random(1, #especies.genomes)]
		g2 = especies.genomes[math.random(1, #especies.genomes)]
		child = crossover(g1, g2)
	else
		g = especies.genomes[math.random(1, #especies.genomes)]
		child = copiaGenome(g)
	end

	mutate(child)

	return child
end

function removeStaleSpecies()
	local survived = {}

	for s = 1,#pool.especies do
		local especies = pool.especies[s]

		table.sort(especies.genomes, function (a,b)
			return (a.distancia > b.distancia)
		end)

		if especies.genomes[1].distancia > especies.melhor_distancia then
			especies.melhor_distancia = especies.genomes[1].distancia
			especies.velha = 0
		else
			especies.velha = especies.velha + 1
		end
		if especies.velha < Especies_Salvas or especies.melhor_distancia >= pool.maxima_distancia then
			table.insert(survived, especies)
		end
	end

	pool.especies = survived
end

function removeWeakSpecies() ---compara as novas especieis com as sobreviventes da iteração anterior para ja exluir as fracas
	local survived = {}

	local sum = Distancia_total_media()
	for s = 1,#pool.especies do
		local especies = pool.especies[s]
		breed = math.floor(especies.Media_distancia / sum * Populacao)
		if breed >= 1 then
			table.insert(survived, especies)
		end
	end

	pool.especies = survived
end


function addToSpecies(child)
	local foundSpecies = false
	for s=1,#pool.especies do
		local especies = pool.especies[s]
		if not foundSpecies and sameSpecies(child, especies.genomes[1]) then
			table.insert(especies.genomes, child)
			foundSpecies = true
		end
	end

	if not foundSpecies then
		local childSpecies = novas_especies()
		table.insert(childSpecies.genomes, child)
		table.insert(pool.especies, childSpecies)
	end
end

function newGeneration() 
	cullSpecies(false) 
	rankGlobally()
	removeStaleSpecies()
	rankGlobally()
	for s = 1,#pool.especies do
		local especies = pool.especies[s]
		Calculo_media_distancia(especies)
	end
	removeWeakSpecies()
	local sum = Distancia_total_media()
	local children = {}
	for s = 1,#pool.especies do
		local especies = pool.especies[s]
		breed = math.floor(especies.Media_distancia / sum * Populacao) - 1
		for i=1,breed do
			table.insert(children, breedChild(especies))
		end
	end
	cullSpecies(true)
	while #children + #pool.especies < Populacao do
		local especies = pool.especies[math.random(1, #pool.especies)]
		table.insert(children, breedChild(especies))
	end
	for c=1,#children do
		local child = children[c]
		addToSpecies(child)
	end

	pool.geracao = pool.geracao + 1

	escreve_poolFile("backup." .. pool.geracao .. "." .. forms.gettext(saveLoadFile))
end

function initializePool()
	pool = novaPool()

	for i=1,Populacao do
		basico = Genome_basico()
		addToSpecies(basico)
	end

	initializeRun()
end

function clearJoypad()
	controller = {}
	for b = 1,#ButtonNames do
		controller["P1 " .. ButtonNames[b]] = false
	end
	joypad.set(controller)
end

function initializeRun() ---função que poe no jogo o indiviuo para fazer o mario andar
	loadState(randomState())
	rightmost = 0
	initialRightMost = 0
	pool.quadro_atual = 0
	timeout = Tempo
	clearJoypad()

	local especies = pool.especies[pool.especie_atual]
	local genome = especies.genomes[pool.genoma_atual]
	generateNetwork(genome)
	evaluateCurrent()
	initialRightMost = marioX
	clock = os.clock()
	frameCount = emu.framecount()
end

function evaluateCurrent() --- avalia a especie no jogo correndo
	local especies = pool.especies[pool.especie_atual]
	local genome = especies.genomes[pool.genoma_atual]

	inputs = getInputs()
	controller = evaluateNetwork(genome.network, inputs)

	if controller["P1 Left"] and controller["P1 Right"] then
		controller["P1 Left"] = false
		controller["P1 Right"] = false
	end
	if controller["P1 Up"] and controller["P1 Down"] then
		controller["P1 Up"] = false
		controller["P1 Down"] = false
	end

	joypad.set(controller)
end

function nextGenome()
	pool.genoma_atual = pool.genoma_atual + 1
	if pool.genoma_atual > #pool.especies[pool.especie_atual].genomes then
		pool.genoma_atual = 1
		pool.especie_atual = pool.especie_atual+1
		if pool.especie_atual > #pool.especies then
			newGeneration()
			pool.especie_atual = 1
		end
	end
end

function Distancia_mensurada() --mede a distancia percorrida e relaciona ela com o respectivo genoma
	local especies = pool.especies[pool.especie_atual]
	local genome = especies.genomes[pool.genoma_atual]

	return genome.distancia ~= 0
end

function displayGenome(genome) ---mostra o genona no boxradius e sua interação 
	local network = genome.network
	local cells = {}
	local i = 1
	local cell = {}
	for dy=-BoxRadius,BoxRadius do
		for dx=-BoxRadius,BoxRadius do
			cell = {}
			cell.x = 50+5*dx
			cell.y = 70+5*dy
			cell.value = network.neurons[i].value
			cells[i] = cell
			i = i + 1
		end
	end
	local biasCell = {}
	biasCell.x = 80
	biasCell.y = 110
	biasCell.value = network.neurons[Inputs].value
	cells[Inputs] = biasCell

	for o = 1,Outputs do
		cell = {}
		cell.x = 220
		cell.y = 30 + 8 * o
		cell.value = network.neurons[MaxNodes + o].value
		cells[MaxNodes+o] = cell
		local color
		if cell.value > 0 then
			color = 0xFF0000FF
		else
			color = 0xFF000000
		end
		gui.drawText(223, 24+8*o, ButtonNames[o], color, 9)
	end

	for n,neuron in pairs(network.neurons) do
		cell = {}
		if n > Inputs and n <= MaxNodes then
			cell.x = 140
			cell.y = 40
			cell.value = neuron.value
			cells[n] = cell
		end
	end

	for n=1,4 do
		for _,gene in pairs(genome.genes) do
			if gene.enabled then
				local c1 = cells[gene.into]
				local c2 = cells[gene.out]
				if gene.into > Inputs and gene.into <= MaxNodes then
					c1.x = 0.75*c1.x + 0.25*c2.x
					if c1.x >= c2.x then
						c1.x = c1.x - 40
					end
					if c1.x < 90 then
						c1.x = 90
					end

					if c1.x > 220 then
						c1.x = 220
					end
					c1.y = 0.75*c1.y + 0.25*c2.y

				end
				if gene.out > Inputs and gene.out <= MaxNodes then
					c2.x = 0.25*c1.x + 0.75*c2.x
					if c1.x >= c2.x then
						c2.x = c2.x + 40
					end
					if c2.x < 90 then
						c2.x = 90
					end
					if c2.x > 220 then
						c2.x = 220
					end
					c2.y = 0.25*c1.y + 0.75*c2.y
				end
			end
		end
	end

	gui.drawBox(50-BoxRadius*5-3,70-BoxRadius*5-3,50+BoxRadius*5+2,70+BoxRadius*5+2,0xFF000000, 0x80808080)
	for n,cell in pairs(cells) do
		if n > Inputs or cell.value ~= 0 then
			local color = math.floor((cell.value+1)/2*256)
			if color > 255 then color = 255 end
			if color < 0 then color = 0 end
			local opacity = 0xFF000000
			if cell.value == 0 then
				opacity = 0x50000000
			end
			color = opacity + color*0x10000 + color*0x100 + color
			gui.drawBox(cell.x-2,cell.y-2,cell.x+2,cell.y+2,opacity,color)
		end
	end
	for _,gene in pairs(genome.genes) do
		if gene.enabled then
			local c1 = cells[gene.into]
			local c2 = cells[gene.out]
			local opacity = 0xA0000000
			if c1.value == 0 then
				opacity = 0x20000000
			end

			local color = 0x80-math.floor(math.abs(Sigmoidal(gene.weight))*0x80)
			if gene.weight > 0 then
				color = opacity + 0x8000 + 0x10000*color
			else
				color = opacity + 0x800000 + 0x100*color
			end
			gui.drawLine(c1.x+1, c1.y, c2.x-3, c2.y, color)
		end
	end

	gui.drawBox(49,71,51,78,0x00000000,0x80FF0000)

	if forms.ischecked(Monstrar_Mutacoes_valores) then
		local pos = 100
		for mutation,rate in pairs(genome.mutacoes_valores) do
			gui.drawText(100, pos, mutation .. ": " .. rate, 0xFF000000, 10)
			pos = pos + 8
		end
	end
end

function getStatePath(name)
	return Diretorio_Armazenamento .. "\\" .. name .. ".state"
end

-- Load/Store states
function loadState(name)
	savestate.load(getStatePath(name))
end

function saveState(name) --- salva o estado
	savestate.save(getStatePath(name))
end

function storeFailure(type) --salva a falha 
	local stateName = "failure\\" .. type .. "\\" .. marioX .. '_' .. marioY
	local statePath = getStatePath(stateName)
	if fileExists(statePath) == false then
		if type == "timedOut" then
			table.insert(estado_pool, stateName)
		end
		saveState(stateName)
	end
end

function storeWIP()
	local stateName = "wip\\" .. marioX .. '_' .. marioY
	local statePath = getStatePath(stateName)
	if fileExists(statePath) == false then
		table.insert(estado_pool, stateName)
		saveState(stateName)
	end
end

function fileExists(name)
	local f= io.open(name, "r")
	if f ~= nil then
		io.close(f)
		return true
	else
		return false
	end
end

function randomState()
	if #estado_pool == 0 then
		console.writeline("Estado não encontrado. Por favor especifique o ultimo treinamento.")
		return "init"
	end
	return estado_pool[math.random(#estado_pool)]
end

function load_estado_pool()
	for file in io.lines(Diretorio_Armazenamento .. "\\training-set.txt", "w") do   ---- training-set.txt
		table.insert(estado_pool, file)
	end
end

function writeestado_pool()
	local estado_poolFile = io.open(Diretorio_Armazenamento .. "\\training-set.txt", "w") ----- usar o arquivo training-set como setar qual arquivo desejamos abrir
----- exemplo, mais de uma fase no mario
	for n, file in pairs(estado_pool) do
		estado_pool_arquivo:write(file .. "\n")
	end
	estado_pool_arquivo:close()
end

function escreve_poolFile(arquivo)
	local file = io.open(Diretorio_Armazenamento .. "\\pool\\" .. arquivo, "w")
	file:write(pool.geracao .. "\n")
	file:write(pool.maxima_distancia .. "\n")
	file:write(#pool.especies .. "\n")
		for n,especies in pairs(pool.especies) do
		file:write(especies.melhor_distancia .. "\n")
		file:write(especies.velha .. "\n")
		file:write(#especies.genomes .. "\n")
		for m,genome in pairs(especies.genomes) do
			file:write(genome.distancia .. "\n")
			file:write(genome.maxneuron .. "\n")
			for mutation,rate in pairs(genome.mutacoes_valores) do
				file:write(mutation .. "\n")
				file:write(rate .. "\n")
			end
			file:write("done\n")

			file:write(#genome.genes .. "\n")
			for l,gene in pairs(genome.genes) do
				file:write(gene.into .. " ")
				file:write(gene.out .. " ")
				file:write(gene.weight .. " ")
				file:write(gene.novo .. " ")
				if(gene.enabled) then
					file:write("1\n")
				else
					file:write("0\n")
				end
			end
		end
		end
		file:close()
end

function savePool()
	local arquivo = forms.gettext(saveLoadFile)
	escreve_poolFile(arquivo)
	writeestado_pool()
end

function loadFile(arquivo)
	local file = io.open(Diretorio_Armazenamento .. "\\pool\\" ..arquivo, "r")
	pool = novaPool()
	pool.geracao = file:read("*number")
	pool.maxima_distancia = file:read("*number")
	forms.settext(max_distancia_Rotulo, "Max distancia: " .. math.floor(pool.maxima_distancia))
	local numSpecies = file:read("*number")
	for s=1,numSpecies do
		local especies = novas_especies()
		table.insert(pool.especies, especies)
		especies.melhor_distancia = file:read("*number")
		especies.velha = file:read("*number")
		local numGenomes = file:read("*number")
		for g=1,numGenomes do
			local genome = newGenome()
			table.insert(especies.genomes, genome)
			genome.distancia = file:read("*number")
			genome.maxneuron = file:read("*number")
			local line = file:read("*line")
			while line ~= "done" do
				genome.mutacoes_valores[line] = file:read("*number")
				line = file:read("*line")
			end
			local numGenes = file:read("*number")
			for n=1,numGenes do
				local gene = newGene()
				table.insert(genome.genes, gene)
				local enabled
				gene.into, gene.out, gene.weight, gene.novo, enabled = file:read("*number", "*number", "*number", "*number", "*number")
				if enabled == 0 then
					gene.enabled = false
				else
					gene.enabled = true
				end

			end
		end
	end
	file:close()

	while Distancia_mensurada() do
		nextGenome()
	end
	initializeRun()
	pool.quadro_atual = pool.quadro_atual + 1
end

function loadPool() ----- carrega a geracao 
	local arquivo = forms.gettext(saveLoadFile)
	loadFile(arquivo)
	load_estado_pool()
end

function logToCsv(pool, distancia) ----- coleta os dados de entrada e saida
	local arquivo = Diretorio_Armazenamento .. "\\" .. Diretorio .. ".csv"
	local file = io.open(arquivo, "a")
	file:write(pool.geracao .. "," .. pool.especie_atual .. "," .. pool.genoma_atual .. "," .. distancia .. "," .. (os.clock() - clock) .. "," .. (emu.framecount() - frameCount) .. "\n")
	file:close()
end




function playTop()
	local maxdistancia = 0
	local maxs, maxg
	for s,especies in pairs(pool.especies) do
		for g,genome in pairs(especies.genomes) do
			if genome.distancia > maxdistancia then
				maxdistancia = genome.distancia
				maxs = s
				maxg = g
			end
		end
	end

	pool.especie_atual = maxs
	pool.genoma_atual = maxg
	pool.maxima_distancia = maxdistancia
	forms.settext(maxdistanciaLabel, "Max distancia: " .. math.floor(pool.maxima_distancia))
	initializeRun()
	pool.quadro_atual = pool.quadro_atual + 1
	return
end

function onExit()
	forms.destroy(form)
end

console.writeline("Novo Treino!")

if pool == nil then
	initializePool()
end

escreve_poolFile("temp.pool")

event.onexit(onExit)

form = forms.newform(200, 260, "distancia")
max_distancia_Rotulo = forms.label(form, "Max distancia: " .. math.floor(pool.maxima_distancia), 5, 8)
showNetwork = forms.checkbox(form, "Monstrar Mapa", 5, 30)
Monstrar_Mutacoes_valores = forms.checkbox(form, "Monstrar Rotas", 5, 52)
restartButton = forms.button(form, "Restart", initializePool, 5, 77)
saveButton = forms.button(form, "Save", savePool, 5, 102)
loadButton = forms.button(form, "Load", loadPool, 80, 102)
saveLoadFile = forms.textbox(form, Diretorio .. ".pool", 170, 25, nil, 5, 148)
saveLoadLabel = forms.label(form, "Save/Load:", 5, 129)
playTopButton = forms.button(form, "Play Top", playTop, 5, 170)
hideBanner = forms.checkbox(form, "Hide Banner", 5, 190)

load_estado_pool()
loadState(randomState())

while true do
	local backgroundColor = 0xD0FFFFFF
	if not forms.ischecked(hideBanner) then
		gui.drawBox(0, 0, 300, 26, backgroundColor, backgroundColor)
	end

	local especies = pool.especies[pool.especie_atual]
	local genome = especies.genomes[pool.genoma_atual]

	if forms.ischecked(showNetwork) then
		displayGenome(genome)
	end

	if pool.quadro_atual%5 == 0 then
		evaluateCurrent()
	end

	if pool.quadro_atual % 250 == 0 then
		storeWIP()
	end

	joypad.set(controller)

	getPositions()
	if marioX > rightmost then ---- evita que o mario corre para direita
		rightmost = marioX
		timeout = Tempo
	end

	timeout = timeout - 1

	local timeoutBonus = pool.quadro_atual / 4
	local isDead = isDead()
	local hasTimedOut = timeout + timeoutBonus <= 0
	if isDead or hasTimedOut then
		local distancia = Funcao_distancia()
		if gameinfo.getromname() == "Super Mario World (USA)" and rightmost > 4816 then
			distancia = distancia + 1000
		elseif gameinfo.getromname() == "Super Mario Bros." and rightmost > 3186 then
			distancia = distancia + 1000
		end
		if distancia == 0 then
			distancia = -1
		end
		genome.distancia = distancia

		if distancia > pool.maxima_distancia then
			pool.maxima_distancia = distancia
			forms.settext(max_distancia_Rotulo, "Max distancia: " .. math.floor(pool.maxima_distancia))
			escreve_poolFile("backup." .. pool.geracao .. "." .. forms.gettext(saveLoadFile))
		end

		logToCsv(pool, distancia) ---- grava os dados em um planilha
		console.writeline("Gene" .. pool.geracao .. " especie " .. pool.especie_atual .. "genoma " .. pool.genoma_atual .. " distancia: " .. distancia)
		pool.especie_atual = 1
		pool.genoma_atual = 1
		while Distancia_mensurada() do
			nextGenome()
		end
		local failureType = isDead and "died" or "timedOut"
		storeFailure(failureType)
		initializeRun()
	end

	local measured = 0
	local total = 0
	for _,especies in pairs(pool.especies) do
		for _,genome in pairs(especies.genomes) do
			total = total + 1
			if genome.distancia ~= 0 then
				measured = measured + 1
			end
		end
	end
	if not forms.ischecked(hideBanner) then
		gui.drawText(0, 0, "Gen " .. pool.geracao .. " especies " .. pool.especie_atual .. " genome " .. pool.genoma_atual .. " (" .. math.floor(measured/total*100) .. "%)" .. " T: " .. timeout + timeoutBonus, 0xFF000000, 10, "Verdana")
		gui.drawText(0, 12, "distancia: " .. string.format("%06i", math.floor(Funcao_distancia())), 0xFF000000, 10, "Verdana")
		gui.drawText(100, 12, "Max distancia: " .. string.format("%06i", math.floor(pool.maxima_distancia)), 0xFF000000, 10, "Verdana")
	end

	pool.quadro_atual = pool.quadro_atual + 1

	emu.frameadvance()
end
